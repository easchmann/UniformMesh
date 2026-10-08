#!/usr/bin/env python3
"""Plots and metrics for a flashx/eval/run.sh result directory.

For each problem the final checkpoints of NewImpl and plain Spark are compared with
the exact solution at the checkpoint's own time (runs stop at nend/tmax, not exactly
at the nominal end time, so the time is read from each checkpoint):

  vortex  isentropic vortex advected by (u_ambient, v_ambient) on the periodic domain;
          the exact solution is the initial vortex moved by (u_ambient*t, v_ambient*t)
  sod     two Riemann problems on the periodic strip, at sim_posn (left|right) and at
          xmin (right|left), solved exactly; valid until their waves meet

Writes <out>/plots/<problem>.png and <out>/results.json, and prints one INFO line per
result (run.sh appends them to summary.txt). The exact solution is checked against the
first checkpoint (the initial condition) before it is used.

usage: plot.py <run.sh output dir> [--compare vortex=4x4 sod=4x1]
needs: numpy, h5py, matplotlib
"""
import argparse
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------- checkpoints

def _named_values(dset):
    return {row["name"].decode().strip().lower(): row["value"] for row in dset[()]}


class Checkpoint:
    """One Flash-X HDF5 checkpoint of a single-block UG run (2D)."""

    def __init__(self, path):
        self.path = Path(path)
        with h5py.File(path, "r") as f:
            if f["bounding box"].shape[0] != 1:
                raise ValueError(f"{path}: expected one block (UG on 1 rank)")
            self.time = float(_named_values(f["real scalars"])["time"])
            self.real = _named_values(f["real runtime parameters"])
            self.int = _named_values(f["integer runtime parameters"])
            # logical parameters are stored as integers (0/1)
            self.bool = {k: bool(v) for k, v in _named_values(f["logical runtime parameters"]).items()}
            names = [n.decode().strip() for n in f["unknown names"][:, 0]]
            # unknowns are (block, k, j, i); keep the k=0 plane as [j, i]
            self.var = {n: np.array(f[n][0, 0]) for n in names}
            bb = np.array(f["bounding box"][0])
        ny, nx = self.var["dens"].shape
        self.lo, self.hi = bb[:2, 0], bb[:2, 1]
        self.dx = (self.hi - self.lo) / [nx, ny]
        self.x = self.lo[0] + (np.arange(nx) + 0.5) * self.dx[0]
        self.y = self.lo[1] + (np.arange(ny) + 0.5) * self.dx[1]


def checkpoints(run_dir):
    return sorted(Path(run_dir).glob("*_hdf5_chk_*"))


# ---------------------------------------------------------------- exact solutions

def vortex_exact(chk, t):
    """Isentropic vortex at time t, as Simulation_initBlock (IsentropicVortex) builds it:
    sub-point averages of rho, rho*u, rho*v, rho*E over nx_subint x ny_subint points per
    cell, centre moved to (xctr + u_ambient*t, yctr + v_ambient*t), nearest periodic image."""
    r, i = chk.real, chk.int
    gamma, beta = r["gamma"], r["vortex_strength"]
    rho0, p0, u0, v0 = r["rho_ambient"], r["p_ambient"], r["u_ambient"], r["v_ambient"]
    nsx, nsy = i.get("nx_subint", 1), i.get("ny_subint", 1)
    length = chk.hi - chk.lo
    xc = chk.lo[0] + (r["xctr"] + u0 * t - chk.lo[0]) % length[0]
    yc = chk.lo[1] + (r["yctr"] + v0 * t - chk.lo[1]) % length[1]
    # sub-points of every cell: shape (ny, nx, nsy, nsx)
    sx = chk.x[None, :, None, None] - 0.5 * chk.dx[0] + (np.arange(nsx)[None, None, None, :] + 0.5) * chk.dx[0] / nsx
    sy = chk.y[:, None, None, None] - 0.5 * chk.dx[1] + (np.arange(nsy)[None, None, :, None] + 0.5) * chk.dx[1] / nsy
    xp = sx - xc
    yp = sy - yc
    xp -= length[0] * np.round(xp / length[0])
    yp -= length[1] * np.round(yp / length[1])
    amp = beta / (2 * math.pi) * np.exp(0.5 * (1 - xp**2 - yp**2))
    tstar = 1 - (gamma - 1) / gamma * 0.5 * amp**2          # T / T_ambient
    rho = rho0 * tstar ** (1 / (gamma - 1))
    p = p0 * tstar ** (gamma / (gamma - 1))
    u = u0 - yp * amp
    v = v0 + xp * amp
    rhoE = p / (gamma - 1) + 0.5 * rho * (u**2 + v**2)
    mean = lambda a: a.mean(axis=(2, 3))
    rho_m = mean(rho)
    u_m, v_m = mean(rho * u) / rho_m, mean(rho * v) / rho_m
    eint = mean(rhoE) / rho_m - 0.5 * (u_m**2 + v_m**2)
    return {"dens": rho_m, "velx": u_m, "vely": v_m, "pres": (gamma - 1) * rho_m * eint}


def riemann_exact(left, right, gamma, xi):
    """Exact solution of the 1D Riemann problem (Toro, ch. 4) at xi = x/t.
    left/right = (rho, u, p); returns rho, u, p arrays and the largest wave speed |s|."""
    rl, ul, pl = left
    rr, ur, pr = right
    cl, cr = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)
    g1, g2 = (gamma - 1) / (2 * gamma), (gamma + 1) / (2 * gamma)

    def f(p, rk, pk, ck):  # pressure function of one side and its derivative
        if p > pk:  # shock
            a, b = 2 / ((gamma + 1) * rk), (gamma - 1) / (gamma + 1) * pk
            s = math.sqrt(a / (p + b))
            return (p - pk) * s, s * (1 - 0.5 * (p - pk) / (b + p))
        # rarefaction
        return 2 * ck / (gamma - 1) * ((p / pk) ** g1 - 1), (p / pk) ** (-g2) / (rk * ck)

    p = max(1e-12, 0.5 * (pl + pr))
    for _ in range(100):
        fl, dl = f(p, rl, pl, cl)
        fr, dr = f(p, rr, pr, cr)
        dp = (fl + fr + ur - ul) / (dl + dr)
        p = max(1e-12, p - dp)
        if abs(dp) < 1e-14 * p:
            break
    us = 0.5 * (ul + ur) + 0.5 * (f(p, rr, pr, cr)[0] - f(p, rl, pl, cl)[0])

    def speed(k):  # outermost wave speed of the left (k=-1) or right (k=+1) wave
        uk, pk, ck = (ul, pl, cl) if k < 0 else (ur, pr, cr)
        if p > pk:
            return uk + k * ck * math.sqrt(g2 * p / pk + g1)  # shock
        return uk + k * ck                                     # rarefaction head

    def side(k, xi):  # sample left (k=-1) or right (k=+1) of the contact
        rk, uk, pk, ck = (rl, ul, pl, cl) if k < 0 else (rr, ur, pr, cr)
        if p > pk:  # shock
            rs = rk * (p / pk + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * p / pk + 1)
            outside = k * (xi - speed(k)) > 0
            return np.where(outside, rk, rs), np.where(outside, uk, us), np.where(outside, pk, p)
        rs = rk * (p / pk) ** (1 / gamma)  # rarefaction
        cs = ck * (p / pk) ** g1
        head, tail = uk + k * ck, us + k * cs
        fan_u = 2 / (gamma + 1) * (-k * ck + (gamma - 1) / 2 * uk + xi)
        # (evaluated everywhere and masked below; clip so points outside the fan stay finite)
        fan_c = np.maximum(2 / (gamma + 1) * ck - k * (gamma - 1) / (gamma + 1) * (uk - xi), 0.0)
        fan_r = rk * (fan_c / ck) ** (2 / (gamma - 1))
        fan_p = pk * (fan_c / ck) ** (2 * gamma / (gamma - 1))
        outside = k * (xi - head) > 0
        inside_fan = ~outside & (k * (xi - tail) > 0)
        r_ = np.where(outside, rk, np.where(inside_fan, fan_r, rs))
        u_ = np.where(outside, uk, np.where(inside_fan, fan_u, us))
        p_ = np.where(outside, pk, np.where(inside_fan, fan_p, p))
        return r_, u_, p_

    xi = np.asarray(xi, dtype=float)
    lft, rgt = side(-1, xi), side(+1, xi)
    is_left = xi < us
    rho, u, pp = (np.where(is_left, a, b) for a, b in zip(lft, rgt))
    return rho, u, pp, max(abs(speed(-1)), abs(speed(+1)))


def sod_exact(real, lo, hi, x, t):
    """Periodic Sod strip at cell centres x: Riemann problems at sim_posn (left|right)
    and at xmin (right|left). Returns the fields and whether the two wave systems are
    still apart (the solution is only exact until they meet)."""
    if abs(real.get("sim_xangle", 0.0)) > 0:
        raise ValueError("sod_exact assumes an interface normal to x (sim_xangle = 0)")
    gamma = real["gamma"]
    left = (real["sim_rholeft"], real["sim_uleft"], real["sim_pleft"])
    right = (real["sim_rhoright"], real["sim_uright"], real["sim_pright"])
    x0, xmin, length = real["sim_posn"], lo[0], hi[0] - lo[0]
    d1 = x - x0
    d2 = x - xmin
    d2 -= length * np.round(d2 / length)
    use1 = np.abs(d1) <= np.abs(d2)
    if t <= 0:
        rho, u, p = (np.where(x < x0, a, b) for a, b in zip(left, right))
        return {"dens": rho, "velx": u, "pres": p}, True
    *s1, smax = riemann_exact(left, right, gamma, d1 / t)
    *s2, _ = riemann_exact(right, left, gamma, d2 / t)
    rho, u, p = (np.where(use1, a, b) for a, b in zip(s1, s2))
    # the two systems meet once their fronts have covered the gap between the interfaces
    gap = min(abs(x0 - xmin), length - abs(x0 - xmin))
    return {"dens": rho, "velx": u, "pres": p}, 2 * smax * t < gap


# ---------------------------------------------------------------- metrics

def norms(a, b):
    d = np.abs(np.asarray(a) - np.asarray(b))
    return {"L1": float(d.mean()), "Linf": float(d.max())}


def mirror_index(chk):
    """Index of the mirror cell about x = (xmin + sim_posn)/2, or None if cells do not map."""
    axis = 0.5 * (chk.lo[0] + chk.real["sim_posn"])
    nx = chk.x.size
    j = (2 * axis - chk.x - chk.lo[0]) / chk.dx[0] - 0.5
    idx = np.round(j).astype(int)
    if np.max(np.abs(j - idx)) > 1e-6:
        return None
    return idx % nx


# ---------------------------------------------------------------- problems

def do_vortex(out, name, ni_dir, ref_dir):
    res, lines = {}, []
    first = Checkpoint(checkpoints(ni_dir)[0])
    init_err = norms(first.var["dens"], vortex_exact(first, first.time)["dens"])
    res["exact_vs_initial_condition"] = init_err
    runs = {"NewImpl": Checkpoint(checkpoints(ni_dir)[-1]), "Spark": Checkpoint(checkpoints(ref_dir)[-1])}
    exact = {}
    for impl, chk in runs.items():
        ex = vortex_exact(chk, chk.time)
        exact[impl] = ex
        res[impl] = {"time": chk.time, "file": str(chk.path),
                     "error_vs_exact": {v: norms(chk.var[v], ex[v]) for v in ("dens", "velx", "vely", "pres")}}
    a, b = runs["NewImpl"], runs["Spark"]
    res["NewImpl_vs_Spark"] = {v: norms(a.var[v], b.var[v]) for v in ("dens", "velx", "vely", "pres")}
    lines.append(f"{name}: exact solution vs stored initial condition: dens Linf {init_err['Linf']:.2e} "
                 f"(check of the exact-solution formula)")
    for impl, chk in runs.items():
        e = res[impl]["error_vs_exact"]["dens"]
        lines.append(f"{name}: {impl} at t={chk.time:.5g}: dens error vs exact L1 {e['L1']:.3e}, Linf {e['Linf']:.3e}")
    # NewImpl and Spark must agree to a small fraction of the discretisation error; with dt
    # from the CFL condition each code picks its own steps, so roundoff-level agreement is
    # not expected (the fixed-dt variant is checked strictly by run.sh)
    diff, err = res["NewImpl_vs_Spark"]["dens"], res["Spark"]["error_vs_exact"]["dens"]
    ratio = {k: diff[k] / err[k] if err[k] > 0 else float("inf") for k in ("L1", "Linf")}
    res["NewImpl_vs_Spark_relative_to_error"] = ratio
    tol = REL_TOL.get(name, REL_TOL_DEFAULT)
    ok = all(r <= tol for r in ratio.values())
    lines.append(f"{'PASS' if ok else 'FAIL'}  {name}: |NewImpl - Spark| / |Spark - exact| (dens) = "
                 f"{ratio['L1']:.2e} (L1), {ratio['Linf']:.2e} (Linf), limit {tol:g}")

    ex = exact["NewImpl"]
    ext = [a.lo[0], a.hi[0], a.lo[1], a.hi[1]]
    fig, ax = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    panels = [
        (ax[0, 0], a.var["dens"], f"NewImpl density, t={a.time:.4g}", "viridis"),
        (ax[0, 1], ex["dens"], "exact density", "viridis"),
        (ax[0, 2], a.var["dens"] - ex["dens"], "NewImpl − exact", "RdBu_r"),
        (ax[1, 0], b.var["dens"], f"Spark density, t={b.time:.4g}", "viridis"),
        (ax[1, 1], b.var["dens"] - exact["Spark"]["dens"], "Spark − exact", "RdBu_r"),
        (ax[1, 2], a.var["dens"] - b.var["dens"], "NewImpl − Spark", "RdBu_r"),
    ]
    for axis, data, title, cmap in panels:
        kw = {}
        if cmap == "RdBu_r":
            m = float(np.max(np.abs(data))) or 1e-300
            kw = {"vmin": -m, "vmax": m}
        im = axis.imshow(data, origin="lower", extent=ext, cmap=cmap, **kw)
        axis.set_title(title)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        fig.colorbar(im, ax=axis, shrink=0.85)
    fig.suptitle(f"{name}: isentropic vortex, density after one period")
    fig.savefig(out / f"{name}.png", dpi=110)
    plt.close(fig)
    return res, lines


def do_sod(out, name, ni_dir, ref_dir):
    res, lines = {}, []
    first = Checkpoint(checkpoints(ni_dir)[0])
    ex0, _ = sod_exact(first.real, first.lo, first.hi, first.x, first.time)
    init_err = norms(first.var["dens"].mean(axis=0), ex0["dens"])
    res["exact_vs_initial_condition"] = init_err
    lines.append(f"{name}: exact solution vs stored initial condition: dens Linf {init_err['Linf']:.2e} "
                 f"(check of the exact-solution setup)")
    runs = {"NewImpl": Checkpoint(checkpoints(ni_dir)[-1]), "Spark": Checkpoint(checkpoints(ref_dir)[-1])}
    prof, exact = {}, {}
    for impl, chk in runs.items():
        ex, valid = sod_exact(chk.real, chk.lo, chk.hi, chk.x, chk.time)
        exact[impl] = ex
        p = {v: chk.var[v].mean(axis=0) for v in ("dens", "velx", "pres")}   # uniform in y
        prof[impl] = p
        y_dev = max(float(np.max(np.abs(chk.var[v] - chk.var[v][0]))) for v in ("dens", "velx", "pres"))
        entry = {"time": chk.time, "file": str(chk.path), "exact_valid": bool(valid),
                 "y_uniformity_max_dev": y_dev,
                 "error_vs_exact": {v: norms(p[v], ex[v]) for v in p},
                 "sum_velx": float(chk.var["velx"].sum())}
        # mirror symmetry about (xmin + sim_posn)/2: u odd, rho even
        mi = mirror_index(chk)
        if mi is not None:
            entry["symmetry"] = {"u_odd_max": float(np.max(np.abs(p["velx"] + p["velx"][mi]))),
                                 "rho_even_max": float(np.max(np.abs(p["dens"] - p["dens"][mi])))}
        res[impl] = entry
        e = entry["error_vs_exact"]["dens"]
        sym = (f"max|u(x)+u(x')| {entry['symmetry']['u_odd_max']:.2e}, "
               f"max|rho(x)-rho(x')| {entry['symmetry']['rho_even_max']:.2e}") if mi is not None else "n/a"
        lines.append(f"{name}: {impl} at t={chk.time:.5g}: dens error vs exact L1 {e['L1']:.3e}, Linf {e['Linf']:.3e}; "
                     f"symmetry {sym}; rows identical in y: {'yes' if y_dev == 0 else f'no (max dev {y_dev:.1e})'}")
        if not valid:
            lines.append(f"{name}: warning: at t={chk.time:.5g} the two periodic wave systems may interact; exact solution not valid")
    res["NewImpl_vs_Spark"] = {v: norms(prof["NewImpl"][v], prof["Spark"][v]) for v in ("dens", "velx", "pres")}

    a = runs["NewImpl"]
    xs = np.linspace(a.lo[0], a.hi[0], 4001)   # exact solution as a curve
    ex_fine, _ = sod_exact(a.real, a.lo, a.hi, xs, a.time)
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    for axis, v, label in [(ax[0, 0], "dens", "density"), (ax[0, 1], "velx", "x-velocity"), (ax[1, 0], "pres", "pressure")]:
        axis.plot(xs, ex_fine[v], "k-", lw=1, label=f"exact, t={a.time:.4g}")
        axis.plot(a.x, prof["NewImpl"][v], "o", ms=3.5, mfc="none", label="NewImpl")
        axis.plot(runs["Spark"].x, prof["Spark"][v], "x", ms=3.5, label="Spark")
        axis.set_title(label)
        axis.set_xlabel("x")
        axis.legend(fontsize=8)
    mi = mirror_index(a)
    axis = ax[1, 1]
    if mi is not None:
        for impl, mk in (("NewImpl", "o-"), ("Spark", "x-")):
            u = prof[impl]["velx"]
            axis.plot(a.x, u + u[mi], mk, ms=3, lw=0.8, label=impl)
        axis.set_title("mirror symmetry: u(x) + u(x′)  (0 if symmetric)")
        axis.set_xlabel("x")
        axis.legend(fontsize=8)
    fig.suptitle(f"{name}: Sod (periodic, 2D strip averaged over y), use_hybridRiemann = {a.bool.get('use_hybridriemann', False)}")
    fig.savefig(out / f"{name}.png", dpi=110)
    plt.close(fig)
    return res, lines


# analysis by problem name prefix: vortex*, sod* (e.g. vortex_fixeddt, sod_hybrid)
KINDS = {"vortex": do_vortex, "sod": do_sod}
REL_TOL = {}            # per problem, from --rel-tol name=value
REL_TOL_DEFAULT = 0.05


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path, help="run.sh output directory")
    ap.add_argument("--compare", nargs="*", default=[], metavar="PROBLEM=SPLIT",
                    help="NewImpl block split to use (default 1x1; all splits are bitwise identical)")
    ap.add_argument("--rel-tol", nargs="*", default=[], metavar="PROBLEM=TOL",
                    help="vortex problems: max |NewImpl - Spark| as a fraction of the error vs exact "
                         f"(default {REL_TOL_DEFAULT})")
    args = ap.parse_args()
    REL_TOL.update({k: float(v) for k, v in (t.split("=", 1) for t in args.rel_tol)})
    split = dict(s.split("=", 1) for s in args.compare)

    plots = args.out / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    results, ok = {}, True
    # every problem directory run.sh wrote (one with a ref/ run), in a stable order
    for name in sorted(d.name for d in args.out.iterdir() if (d / "ref").is_dir()):
        func = next((f for kind, f in KINDS.items() if name.startswith(kind)), None)
        if func is None:
            print(f"INFO  {name}: no analysis for this problem, skipped")
            continue
        ni, ref = args.out / name / f"ni_{split.get(name, '1x1')}", args.out / name / "ref"
        if not (checkpoints(ni) and checkpoints(ref)):
            print(f"INFO  {name}: no checkpoints in {ni} or {ref}, skipped")
            continue
        try:
            res, lines = func(plots, name, ni, ref)
        except Exception as exc:  # report and continue with the next problem
            print(f"INFO  {name}: plot.py failed: {exc}")
            ok = False
            continue
        res["plot"] = str(plots / f"{name}.png")
        results[name] = res
        for line in lines:   # checks carry their own PASS/FAIL, everything else is INFO
            print(line if line.startswith(("PASS  ", "FAIL  ")) else f"INFO  {line}")
        print(f"INFO  {name}: plot {plots / f'{name}.png'}")
    (args.out / "results.json").write_text(json.dumps(results, indent=2))
    print(f"INFO  metrics in {args.out / 'results.json'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
