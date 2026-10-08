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
                                       # (ONLY="vortex_fixeddt" for selected problems)
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

## Problems (`vortex.par`, `vortex_fixeddt.par`, `sod.par`, `sod_hybrid.par`)

Both are Flash-X's own simulation units (`Simulation/SimulationMain/IsentropicVortex`,
`.../Sod`); the parameter files here adapt them to the mesh (uniform blocks, periodic).

| | Isentropic vortex | Sod | Sod, hybrid solver |
|---|---|---|---|
| Grid | 256 x 256, domain [0,10]^2 | 256 x 256, domain [0,1]^2 (uniform in y) | same as Sod |
| Boundaries | periodic | periodic (two mirror-image wave systems: at x = 0.5 and at x = 0) | same |
| End | t = 10 exactly (one period) | t = 0.025 | same |
| dt | set by the CFL condition (`dtmax` large, `nend` large) | same | same |
| CFL | 0.5 | 0.5 | 0.5 |
| Riemann solver | HLLC (no shocks) | HLLC everywhere (`use_hybridRiemann = .false.`, Spark's default) | HLLE in shock-flagged cells (`use_hybridRiemann = .true.`) |
| What it tests | smooth flow | shocks: shock detection, flattening, viscosity | the hybrid HLLE path in addition |

`vortex_fixeddt` is the vortex with a fixed dt = 0.005 (below the CFL limit of ~0.0067), so
NewImpl and Spark take identical steps; it runs on the vortex executables. With dt from the
CFL condition each code picks dt from its own solution, and roundoff-level differences
between the two codes change the step sequence, so only the fixed-dt variant can be
compared to roundoff. `sod_hybrid` runs on the Sod executables; only the parameter file differs. NewImpl takes
`use_hybridRiemann`, `use_flattening`, `cvisc`, `smlrho` and `smallp` from Spark's runtime
parameters (as NewImpl's own `hy_rk_getFaceFlux_wrapper` does), so both codes always run
with the same solver settings.

Differences to Flash-X's standard Sod (`Sod/flash.par`): a diagonal interface
(`sim_xangle = sim_yangle = 45`) on AMR with outflow boundaries to t = 0.2. Here the interface
is normal to x on a uniform grid with periodic boundaries (the mesh supports only those),
which adds the mirrored problem at x = 0 and limits the run to the time before the two
wave systems meet.

## Checks (`run.sh`)

1. **Block-split invariance**: vortex with blocks 2x2, 4x4, 8x8, 32x32, 4x1, 1x16,
   vortex_fixeddt with 4x4, 32x32, Sod with 4x1, 8x2, 16x16, 32x4 and Sod-hybrid with 4x1,
   16x16; every final checkpoint must be bitwise
   identical to the 1x1 run (`sfocu` SUCCESS). Tests halo exchange, flux synchronisation
   and the neighbour lookup. Every block has at least 8 cells per axis (`MIN_BLOCK_CELLS`
   in `config.sh`; `run.sh` refuses smaller splits).
2. **NewImpl vs plain Spark**:
   - `vortex_fixeddt` (identical steps) must agree to roundoff: max `sfocu` mag error
     <= 1e-12 (`VORTEX_TOL`);
   - `vortex` (dt from CFL) must agree to a small fraction of the discretisation error:
     |NewImpl − Spark| <= 5 % of |Spark − exact| in L1 and Linf (`REL_TOL`, checked by `plot.py`);
   - the Sod problems are reported only (INFO): NewImpl and Spark differ near shocks, and
     with dt from CFL the two runs end at different times, so they are judged against the
     exact solution; the summary lists the error norms and the sum of `velx`, which is 0
     for a solution that keeps the problem's mirror symmetry about x = 0.25.
3. Both take the same number of steps.

4. **Exact solution and plots** (`plot.py`; information, plus the `vortex` REL_TOL check): the final checkpoints of
   NewImpl and Spark are compared with the exact solution **at the checkpoint's own time**
   (in general a run need not end exactly at the nominal time, e.g. when it stops at
   `nend`, so the exact solution is evaluated at the stored time):
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
├── plots/sod_hybrid.png, plots/vortex_fixeddt.png  the same for the variants
├── plot.log           plot.py output (also appended to summary.txt)
└── vortex/, vortex_fixeddt/, sod/, sod_hybrid/  ni_<split>/ and ref/ run directories (flash.par, run.out, *.log,
                       *_hdf5_chk_*, *_hdf5_plt_cnt_*), sfocu_*.txt reports
```

`plot.py` needs numpy, h5py and matplotlib (on riken: `python3 -m pip install --user "h5py<3.13"`);
without them `run.sh` skips it, and it can be run later: `flashx/eval/plot.py <eval_runs/date>`.
Checkpoints open in VisIt/yt, or with `h5dump`. To look at a result on the laptop:

```sh
rsync -az --exclude '*_hdf5_*' riken:Flash-X-MOL/eval_runs/<date>/ eval_runs/<date>/   # without the HDF5 files
```

Results (Flash-X `44188e65b`, UniformMesh `c696ce3` + the solver-settings fix; obtained
with the earlier settings: vortex 64x64 and Sod 128x16, fixed or capped dt, CFL 0.8, splits
down to 4 cells per block):
- all block splits bitwise identical (vortex, Sod, Sod-hybrid);
- vortex: NewImpl = Spark to ~1e-14; error vs exact at t = 9.979: dens L1 3.3e-4,
  Linf 7.6e-3 (identical for both);
- Sod at t = 0.0255, dens error vs exact (L1 / Linf) and mirror symmetry max |u(x)+u(x')|:

  | | NewImpl | Spark |
  |---|---|---|
  | `sod` (HLLC everywhere) | 6.05e-3 / 7.24e-2, symmetry 1.5e-2 | 6.41e-3 / 7.54e-2, symmetry 1e-14 |
  | `sod_hybrid` (HLLE in shocks) | 6.08e-3 / 7.23e-2, symmetry 1.5e-2 | 6.44e-3 / 7.59e-2, symmetry 1e-14 |

  NewImpl breaks the mirror symmetry at the waves with and without the hybrid solver;
  Spark keeps it in both. The standalone NewImpl driver, without the UniformMesh, breaks it too.

## Not covered

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
