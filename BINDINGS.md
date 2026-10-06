# Fortran bindings

How Flash-X Fortran code (e.g. Spark's `Hydro.F90`) uses the C mesh.
This is a thin binding: a Fortran module that mirrors the C API. It does not implement Flash-X's `Grid_interface`; the Spark driver calls the module directly instead of `Grid_*`.

Files: `mesh_bind.c` (C side), `mesh_f.F90` (module `mesh_f`), `test_mesh_f.F90` (tests).

---

## Design

```
Flash-X (Fortran)        use mesh_f
                             |   1-based IDs, global indices, real(mesh_rk) pointers
mesh_f.F90  ---------------  |   bind(C) interfaces, layout detection, index shift
                             |   plain C ABI: ints, doubles, pointers
mesh_bind.c / mesh.c  -----  |   uses the MESH_* macros
mesh_config.h               glue: layout, dimension, halo, type
```

Rules:
1. **Fortran never sees the `Mesh` struct or the macros.** The handle is an opaque `type(c_ptr)` inside `mesh_t`. All index arithmetic (`MESH_GIDX`, `MESH_BCOORD`, `MESH_IDX`, `MESH_X`) stays in C, so the glue stays the single source of truth.
2. **Zero copy.** `mesh_f_data_ptr` returns a Fortran pointer that aliases the C storage of one block. Kernels read and write the mesh directly.
3. **Flash-X conventions on the Fortran side**, so Spark code works with its usual indexing.
4. Calls through the binding are not on the hot path (a few per block). Hot loops only touch the pointer.

---

## Conventions

| | C | Fortran |
|---|---|---|
| block | `b = 0 .. NBLOCK-1` | `blockID = 1 .. nblocks` |
| variable | `v = 0 .. NVAR-1` | `var = 1 .. nvars` (e.g. `DENS_VAR`) |
| axis | `a = 0, 1, 2` | `1, 2, 3` (`IAXIS..KAXIS`) |
| cell index | local, interior `0..N-1` | **global**, interior of the domain `1..N_global` |
| unused axis | `0..0` | `1..1` |
| block limits | — | `limits(LOW:HIGH, MDIM)`, like `blkLimits` / `blkLimitsGC` |

The Fortran cell index is `MESH_GIDX(m,b,a,idx) + 1`. This is the AMReX-style global index. With it the Flash-X coordinate formulas hold unchanged with `x0 = domain low`:
`center = x0 + (i - 0.5)*dx`, `left = x0 + (i - 1)*dx`, `right = x0 + i*dx` (`Grid_getCenterCoords` etc.).

`MESH_LOW = 1`, `MESH_HIGH = 2`, `MESH_LEFT_EDGE = 1`, `MESH_CENTER = 2` and `MESH_RIGHT_EDGE = 3` have the same values as Flash-X's `LOW`, `HIGH`, `LEFT_EDGE`, `CENTER` and `RIGHT_EDGE` (`constants.h`), so the arrays and arguments can be passed as they are.

---

## Data layout

The layout is still decided only by `MESH_IDX` in the glue. `mesh_bind.c` reports the address of a block's first guard cell and the element strides in `v, i, j, k`. `mesh_f_create` accepts the two orders that a contiguous Fortran array can represent:

| C layout | Strides `(v, i, j, k)` | Fortran pointer | Flash-X equivalent |
|---|---|---|---|
| AoS | `1, nvar, nvar*nt1, nvar*nt1*nt2` | `U(1:nvar, iloGC:, jloGC:, kloGC:)` | `UG` (`unk(var,i,j,k)`) |
| SoA | `ncell, 1, nt1, nt1*nt2` | `U(iloGC:, jloGC:, kloGC:, 1:nvar)` | `UGReordered`, `!!Reorder(4)` |

Any other `MESH_IDX` makes `mesh_f_create` fail with `MESH_ERR_LAYOUT`. Both layouts keep each block contiguous (block outermost), which is what makes the per-block pointer possible.

The order of the Fortran array must match the order the Flash-X code was compiled for. Spark marks its arrays with `!!Reorder(4)`, and Flash-X's setup rewrites them for the reordered layout. Check `mesh_f_layout(m)` against it at init (see below).

---

## Types

`mesh_rk` is the kind of the mesh data. It is chosen at compile time from `MESH_TYPE_BYTES` in `mesh_config.h` (`c_double` for 8, `c_float` for 4). `mesh_f_create` also checks the size against the C `sizeof(type_t)` at run time.

Flash-X compiles with `real` promoted to 8 bytes, so `real` and `real(mesh_rk)` agree for a `double` mesh. A `float` mesh would need Flash-X built with 4-byte `real`. Coordinates and deltas are always `real(c_double)`.

---

## API (`use mesh_f`)

| Procedure | Purpose |
|---|---|
| `mesh_f_create(m, nx, ny, nz, bx, by, bz, nvar, ierr)` | global cells and blocks per axis; `ierr = MESH_OK` or `MESH_ERR_CREATE / _TYPE / _LAYOUT` |
| `mesh_f_remove(m)` | free the mesh |
| `mesh_f_ndim()` | `MESH_NDIM` of the C side |
| `mesh_f_nblocks(m)`, `mesh_f_nvars(m)`, `mesh_f_nhalo(m)` | sizes |
| `mesh_f_layout(m)` | `MESH_LAYOUT_AOS` or `MESH_LAYOUT_SOA` |
| `mesh_f_limits(m, blockID, limits [, limitsGC])` | interior / guard-cell limits, global 1-based |
| `mesh_f_data_ptr(m, blockID, U)` | pointer to the block's data incl. halo (zero copy) |
| `mesh_f_set_domain(m, low, high, ierr)` | set the global domain (default `[0,1]` per axis); `dx` follows. Must match the Flash-X domain, because Spark takes `dx` from the mesh |
| `mesh_f_deltas(m, deltas)` | `dx, dy, dz` |
| `mesh_f_domain(m, low, high)` | global domain corners (`low` is Spark's `hy_globalLBnd`) |
| `mesh_f_cell_coords(m, blockID, axis, edge, coords)` | centers / edges of all cells of the block incl. halo, like `Grid_getCellCoords` |
| `mesh_f_fill_halo(m)` | periodic halo fill, all variables |
| `mesh_f_fill_halo_vars(m, var_first, nvars, ierr)` | halo fill of a variable range (1-based) |
| `mesh_f_write_dump(m, filename, names, step, time, ierr)` | binary dump of all blocks incl. halo (format in `mesh_f.F90`), plotted by `flashx/plot_mesh.py` |
| `mesh_f_c_ptr(m)` | raw C handle for calling other C functions |

---

## Use in Spark's `Hydro.F90`

The block loop replaces the tile iterator. The scratch mapping, `Hydro_prepBlock` and `Hydro_advance` stay as they are.

```fortran
use mesh_f
type(mesh_t), save :: hy_mesh          ! created once at init

! init: create and check consistency with the Flash-X build
call mesh_f_create(hy_mesh, nx, ny, nz, bx, by, bz, NUNK_VARS, ierr)
if (ierr /= MESH_OK)                 call Driver_abort("mesh_f_create failed")
if (mesh_f_ndim()  /= NDIM)          call Driver_abort("NDIM mismatch")
if (mesh_f_nhalo(hy_mesh) /= NGUARD) call Driver_abort("NGUARD /= MESH_NHALO")
if (mesh_f_layout(hy_mesh) /= MESH_LAYOUT_AOS) call Driver_abort("layout /= unk(var,i,j,k)")  ! SOA for a reordered build
call mesh_f_domain(hy_mesh, hy_globalLBnd, high)

! each step (telescoping): one halo fill, then all blocks
call mesh_f_fill_halo(hy_mesh)       ! replaces @M hy_globalFillGuardCells
call mesh_f_deltas(hy_mesh, deltas)
do blockID = 1, mesh_f_nblocks(hy_mesh)
   call mesh_f_limits(hy_mesh, blockID, blkLimits, blkLimitsGC)
   call mesh_f_data_ptr(hy_mesh, blockID, Uin)
   lo(:)   = blkLimits(LOW, :)
   loGC(:) = blkLimitsGC(LOW, :)
   ! ... @M hy_map_scr_ptrs, Hydro_prepBlock(Uin, ...), Hydro_advance(stage, Uin, ...) unchanged
   nullify(Uin)
end do
```

Notes:
- `Grid_fillGuardCells(..., doEos=.TRUE.)` needs no EOS call here. Periodic guard cells are copies of interior cells that are already EOS-consistent.
- Flux correction (`Grid_putFluxData`, `hy_rk_correctFluxes`) is only for AMR, so keep `hy_fluxCorrect = .false.`.
- `MESH_NHALO` must be `MAXSTAGE * NSTENCIL` (Spark's `GUARDCELLS`), e.g. `macros/halo_8.ini` for RK2 with WENO/MP5 (see MACROS.md).
- `Uin` must be a plain `pointer` (as in `Hydro.F90`). Spark's kernels take it as an assumed-shape dummy `dimension(1:, loGC(1):, loGC(2):, loGC(3):)`, which matches the pointer bounds.

---

## Build and test

```
make check DIM=2d LAYOUT=soa HALO=1   # one variant: C tests and Fortran tests
make check-all                        # 2d/3d x soa/aos x halo 1/4/8
```

`test_mesh_f` is built with `-fcheck=all` (bounds checking) and checks:
- the blocks' limits tile the global index space exactly once;
- pointer bounds; a value written through the Fortran pointer is read back by C `mesh_get` at the matching local index, for every cell and variable;
- `mesh_f_fill_halo` and `mesh_f_fill_halo_vars` through the binding: the halo holds the periodic partner, and the halos of unfilled variables stay untouched;
- coordinates against the Flash-X formulas above;
- a Spark-shaped kernel (assumed-shape dummy with `loGC` lower bounds, loops over the grown interior, `±1` neighbour reads) against exact expected values.

The test is compiled with `-DMESH_SOA` for the SoA layout, to pick the array order. Flash-X does the same with `!!Reorder`.


---

## Not covered (yet)

- `Grid_interface` itself (`Grid_getTileIterator`, `Grid_tile_t`, `Grid_fillGuardCells`, ...). The driver calls `mesh_f` directly.
- Scratch arrays (`hy_starState`, fluxes, ...) still come from Spark's own allocation.
- MPI, tiling within a block, AMR, GPU data movement.
- Choosing scattered variables for the halo fill (Spark's `hy_gcMask`). Only a contiguous range is supported.
