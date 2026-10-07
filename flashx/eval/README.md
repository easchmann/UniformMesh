# Evaluating UniformMeshNewImpl

Scripts that reproduce the builds and the evaluation of NewImpl on the UniformMesh
(the `UniformMeshNewImpl` Spark sub-unit) on the RIKEN cluster, from scratch:

```sh
git clone -b MOL git@github.com:RIKEN-RCCS/Flash-X.git ~/Flash-X-MOL     # or an existing checkout
git -C ~/Flash-X-MOL checkout 44188e65bc44fac1982b803aa38cc6d5acfa244e
git clone https://github.com/easchmann/UniformMesh ~/UniformMesh

~/UniformMesh/flashx/eval/prepare.sh   # once: HDF5, site file, Spark/Makefile fix, sfocu
~/UniformMesh/flashx/eval/build.sh     # 4 executables
~/UniformMesh/flashx/eval/run.sh       # runs + comparisons -> eval_runs/<date>/summary.txt
```

All paths and versions are in `config.sh` and can be overridden from the environment
(`FLASHX`, `MPI_PATH`, `HDF5_PAR`, `HDF5_SER`, `JOBS`, `OUT`, ...).

## What gets built

| Object directory (in `$FLASHX`) | Problem | Hydro |
|---|---|---|
| `eval_vortex_ni` | IsentropicVortex | NewImpl on the UniformMesh |
| `eval_vortex_ref` | IsentropicVortex | plain Spark |
| `eval_sod_ni` | Sod | NewImpl on the UniformMesh |
| `eval_sod_ref` | Sod | plain Spark |

All are `./setup <problem> -auto -2d +spark +nofbs -site=supercomp` (UG Grid, 1 MPI rank)
with HDF5 IO. The NewImpl builds add `--with-unit=physics/Hydro/HydroMain/Spark/UniformMeshNewImpl`;
`build.sh` reinstalls that sub-unit from this repo (`flashx/install_newimpl.sh`, dim 2d,
halo 4) before building, so they always use the current working copy.

## Problems (`vortex.par`, `sod.par`)

| | Isentropic vortex | Sod |
|---|---|---|
| Grid | 64 x 64, domain [0,10]^2 | 128 x 16, domain [0,1]^2 |
| Boundaries | periodic | periodic (two mirror-image wave systems: at x = 0.5 and at x = 0) |
| End | t = 10 (one period, 400 steps) | t = 0.025 (26 steps) |
| dt | fixed 0.025 (`dtmin = dtmax`) | capped at 1e-3 (`dtmax`) |
| What it tests | smooth flow | shocks: shock detection, hybrid Riemann solver, flattening |

## Checks (`run.sh`)

1. **Block-split invariance**: vortex with blocks 2x2, 4x4, 8x8, 16x16, 4x1, 1x8 and Sod
   with 4x1, 8x2, 16x4, 32x1; every final checkpoint must be bitwise identical to the
   1x1 run (`sfocu` SUCCESS). Tests halo exchange, flux synchronisation and the
   neighbour lookup.
2. **NewImpl vs plain Spark**: vortex must agree to roundoff (max `sfocu` mag error
   <= 1e-12). Sod is reported only (INFO): NewImpl and Spark differ near shocks; the
   summary lists the error norms and the sum of `velx`, which is 0 for a solution that
   keeps the problem's mirror symmetry about x = 0.25.
3. Both take the same number of steps.

4. **Exact solution and plots** (`plot.py`, information only): the final checkpoints of
   NewImpl and Spark are compared with the exact solution **at the checkpoint's own time**
   (runs stop at `nend` or just past `tmax`, e.g. the vortex at t = 9.979, not 10, so
   "last minus first checkpoint" would include the vortex's missing 0.021 of travel):
   - vortex: the initial vortex moved by (u_ambient·t, v_ambient·t), built exactly like
     `Simulation_initBlock` (sub-point averages, nearest periodic image);
   - Sod: exact Riemann solutions (Toro) at `sim_posn` (left|right) and at `xmin`
     (right|left); flagged if the two periodic wave systems may already interact.

   Both are first checked against the stored initial condition (vortex ~4e-16, Sod 0).
   Reported: L1/Linf errors vs exact, NewImpl vs Spark, the Sod mirror symmetry
   (max |u(x)+u(x')|, max |rho(x)-rho(x')|) and whether all rows in y are identical.

## Results and inspection

```
eval_runs/<date>/
├── summary.txt        one PASS/FAIL/INFO line per check (run.sh exits with 1 if a check fails)
├── versions.txt       date, host, commits (+ uncommitted changes), HDF5, compiler, dim/halo
├── results.json       all plot.py metrics, per problem and implementation
├── plots/vortex.png   density NewImpl / exact / Spark, errors vs exact, NewImpl − Spark
├── plots/sod.png      rho, u, p profiles vs exact; mirror symmetry u(x) + u(x')
└── vortex/, sod/      ni_<split>/ and ref/ run directories (flash.par, run.out, *.log,
                       *_hdf5_chk_*, *_hdf5_plt_cnt_*), sfocu_*.txt reports
```

`plot.py` needs numpy, h5py and matplotlib (on riken: `python3 -m pip install --user "h5py<3.13"`);
without them `run.sh` skips it, and it can be run later: `flashx/eval/plot.py <eval_runs/date>`.
Checkpoints open in VisIt/yt, or with `h5dump`. To look at a result on the laptop:

```sh
rsync -az --exclude '*_hdf5_*' riken:Flash-X-MOL/eval_runs/<date>/ eval_runs/<date>/   # without the HDF5 files
```

Result of the first evaluation (Flash-X `44188e65b`, UniformMesh `c696ce3`):
- all block splits bitwise identical;
- vortex: NewImpl = Spark to ~1e-14; error vs exact at t = 9.979: dens L1 3.3e-4,
  Linf 7.6e-3 (identical for both);
- Sod at t = 0.0255: dens error vs exact L1 6.1e-3 (NewImpl) vs 6.4e-3 (Spark), Linf
  7.2e-2 vs 7.5e-2. Spark keeps the mirror symmetry to 1e-14, NewImpl breaks it at the
  waves by up to 1.5e-2 (u) and 1.0e-2 (rho). The standalone NewImpl driver, without the UniformMesh, breaks it too.
## Not covered

- dt is fixed/capped in both problems, so the CFL-limited path (and the dt tolerance in
  `Hydro.F90-mc`) is not exercised; a convergence study with a large `dtmax` would be.
- Order of accuracy (one resolution per problem), RK3 (`Spark/rk3`) builds, 3D, more
  than one MPI rank.
- Restart from a checkpoint: the UG reader of `hdf5/parallel/NoFbs` calls the stub
  `Grid_getBlkIndexLimits`, so restarts do not work; writing checkpoints does.

## Environment notes

- `prepare.sh` changes the Flash-X checkout in two places: `sites/supercomp/Makefile.h`
  (a previous version is backed up) and `source/physics/Hydro/HydroMain/Spark/Makefile`,
  where MOL lists `hy_getFaceFlux.o` and `hy_updateSolution.o` without sources (plain
  Spark does not build on MOL otherwise; the original is kept as `Makefile.orig`).
- HDF5 must be parallel for Flash-X: with `+nofbs` the UG Grid selects
  `IO/IOMain/hdf5/parallel/NoFbs`; the serial UG writer needs a fixed block size.
- Always call `./setup`, not `./bin/setup.py` from the checkout root: the latter does
  not find `setup_shortcuts.txt` and ignores every `+shortcut` without an error.
