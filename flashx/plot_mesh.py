#!/usr/bin/env python3
"""Plot UniformMesh dumps (um_dump_<step>.bin, see mesh_f_write_dump in mesh_f.F90).

One dump:   the field over the domain with the mesh blocks drawn in, and one block
            with its halo (the halo cells are copies of the neighbouring blocks).
Two dumps:  both fields and their difference, e.g.
              - initial vs final state of the isentropic vortex (it returns to its start
                at t = 10, so the difference is the numerical error), or
              - two runs with different block splits at the same step (should be 0).

usage:
  python3 plot_mesh.py um_dump_000000.bin                        # field + mesh blocks
  python3 plot_mesh.py um_dump_000000.bin um_dump_000400.bin     # initial vs final
  python3 plot_mesh.py run4x4/um_dump_000400.bin run8x8/um_dump_000400.bin --labels "4x4" "8x8"
options: --var dens (variable), --block N (block shown with its halo), --k K (slice in 3D),
         -o file.png (default: mesh.png)

Needs numpy and matplotlib.
"""

import argparse
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

# colors (reference data-viz palette, light mode)
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SEQUENTIAL = LinearSegmentedColormap.from_list(   # one hue, light -> dark (blue 100 -> 700)
    "blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIVERGING = LinearSegmentedColormap.from_list(    # blue <-> gray midpoint <-> red
    "blue_red", ["#184f95", "#6da7ec", "#f0efec", "#ec8a89", "#a52a2a"])
BLOCK_LINE = "#8a8984"                             # block boundaries: neutral, visible on light and dark cells
HALO_SHADE = (0.99, 0.99, 0.98, 0.35)              # light veil over the halo cells of the block view


def read_dump(path):
    """Read a dump into a dict: header values, variable names and per-block arrays [var, k, j, i]."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:8] != b"UMDUMP01":
        sys.exit(f"error: {path} is not a UniformMesh dump")
    pos = 8

    def ints(n):
        nonlocal pos
        a = np.frombuffer(raw, dtype=np.int32, count=n, offset=pos)
        pos += 4 * n
        return a

    def doubles(n):
        nonlocal pos
        a = np.frombuffer(raw, dtype=np.float64, count=n, offset=pos)
        pos += 8 * n
        return a

    ndim, nvar, nblock, step, layout, nbytes = ints(6)
    gsize, nb, nh = ints(3), ints(3), ints(3)
    time = doubles(1)[0]
    low, high = doubles(3), doubles(3)
    names = []
    for _ in range(nvar):
        names.append(raw[pos:pos + 8].decode().strip())
        pos += 8

    dtype = np.float64 if nbytes == 8 else np.float32
    blocks = []
    for _ in range(nblock):
        lim = ints(6).reshape(3, 2)                # lim[axis] = (low, high), global 1-based
        nt = lim[:, 1] - lim[:, 0] + 1 + 2 * nh    # cells per axis incl. halo
        count = int(nvar * nt.prod())
        data = np.frombuffer(raw, dtype=dtype, count=count, offset=pos)
        pos += nbytes * count
        if layout == 1:   # AoS, Fortran (nvar, i, j, k)
            data = data.reshape(nt[2], nt[1], nt[0], nvar).transpose(3, 0, 1, 2)
        else:             # SoA, Fortran (i, j, k, nvar)
            data = data.reshape(nvar, nt[2], nt[1], nt[0])
        blocks.append({"lim": lim, "data": data})

    return {"path": path, "ndim": int(ndim), "nvar": int(nvar), "step": int(step),
            "time": float(time), "gsize": gsize, "nb": nb, "nh": nh,
            "low": low, "high": high, "names": names, "blocks": blocks}


def var_index(d, var):
    if var not in d["names"]:
        sys.exit(f"error: variable '{var}' not in {d['path']} (has: {', '.join(d['names'])})")
    return d["names"].index(var)


def global_field(d, v, k):
    """Interior of all blocks assembled into one array [j, i] (slice k for 3D)."""
    g = np.full((d["gsize"][1], d["gsize"][0]), np.nan)
    nh = d["nh"]
    for b in d["blocks"]:
        lo, hi = b["lim"][:, 0], b["lim"][:, 1]
        if not (lo[2] <= k <= hi[2]):
            continue
        kk = k - lo[2] + nh[2]
        interior = b["data"][v, kk, nh[1]:nh[1] + hi[1] - lo[1] + 1, nh[0]:nh[0] + hi[0] - lo[0] + 1]
        g[lo[1] - 1:hi[1], lo[0] - 1:hi[0]] = interior
    return g


def cell_width(d):
    return (d["high"] - d["low"]) / d["gsize"]


def draw_blocks(ax, d, k):
    """Block boundaries over the domain."""
    dx = cell_width(d)
    for b in d["blocks"]:
        lo, hi = b["lim"][:, 0], b["lim"][:, 1]
        if not (lo[2] <= k <= hi[2]):
            continue
        x0, y0 = d["low"][0] + (lo[0] - 1) * dx[0], d["low"][1] + (lo[1] - 1) * dx[1]
        ax.add_patch(Rectangle((x0, y0), (hi[0] - lo[0] + 1) * dx[0], (hi[1] - lo[1] + 1) * dx[1],
                               fill=False, edgecolor=BLOCK_LINE, linewidth=0.9))


def style(ax, title):
    ax.set_title(title, loc="left", color=TEXT_PRIMARY, fontsize=10)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    ax.set_xlabel("x", color=TEXT_SECONDARY)
    ax.set_ylabel("y", color=TEXT_SECONDARY)
    for s in ax.spines.values():
        s.set_visible(False)


def colorbar(fig, im, ax, label):
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(label, color=TEXT_SECONDARY, fontsize=9)
    cb.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    cb.outline.set_visible(False)


def domain_extent(d):
    return [d["low"][0], d["high"][0], d["low"][1], d["high"][1]]


def mesh_summary(d):
    n = d["gsize"] // d["nb"]
    nd = d["ndim"]
    blocks = " x ".join(str(x) for x in d["nb"][:nd])
    cells = " x ".join(str(x) for x in n[:nd])
    return f"{blocks} blocks of {cells} cells, halo {d['nh'][0]}"


def plot_single(d, var, block, k, out):
    v = var_index(d, var)
    g = global_field(d, v, k)
    fig, (ax_f, ax_b) = plt.subplots(1, 2, figsize=(13, 5.8), facecolor=SURFACE,
                                     gridspec_kw={"width_ratios": [1.15, 1]})

    # field over the domain with the block boundaries
    im = ax_f.imshow(g, origin="lower", extent=domain_extent(d), cmap=SEQUENTIAL, interpolation="nearest")
    draw_blocks(ax_f, d, k)
    colorbar(fig, im, ax_f, var)
    style(ax_f, f"{var} at t = {d['time']:.4g} (step {d['step']})\n{mesh_summary(d)}")

    # one block with its halo; the interior is outlined, the halo is veiled
    if block is None:   # default: the block with the smallest value (e.g. the vortex core)
        block = int(np.argmin([np.nanmin(b["data"][v]) for b in d["blocks"]])) + 1
    if not 1 <= block <= len(d["blocks"]):
        sys.exit(f"error: --block must be between 1 and {len(d['blocks'])}")
    b = d["blocks"][block - 1]
    lo, hi, nh = b["lim"][:, 0], b["lim"][:, 1], d["nh"]
    kk = min(max(k, lo[2]), hi[2]) - lo[2] + nh[2]
    dx = cell_width(d)
    x0 = d["low"][0] + (lo[0] - 1 - nh[0]) * dx[0]
    x1 = d["low"][0] + (hi[0] + nh[0]) * dx[0]
    y0 = d["low"][1] + (lo[1] - 1 - nh[1]) * dx[1]
    y1 = d["low"][1] + (hi[1] + nh[1]) * dx[1]
    im_b = ax_b.imshow(b["data"][v, kk], origin="lower", extent=[x0, x1, y0, y1], cmap=SEQUENTIAL,
                       interpolation="nearest", vmin=np.nanmin(g), vmax=np.nanmax(g))
    ix0 = d["low"][0] + (lo[0] - 1) * dx[0]
    iy0 = d["low"][1] + (lo[1] - 1) * dx[1]
    iw, ih = (hi[0] - lo[0] + 1) * dx[0], (hi[1] - lo[1] + 1) * dx[1]
    # veil the halo: four strips around the interior
    for rx, ry, rw, rh in [(x0, y0, x1 - x0, iy0 - y0), (x0, iy0 + ih, x1 - x0, y1 - iy0 - ih),
                           (x0, iy0, ix0 - x0, ih), (ix0 + iw, iy0, x1 - ix0 - iw, ih)]:
        ax_b.add_patch(Rectangle((rx, ry), rw, rh, facecolor=HALO_SHADE, edgecolor="none"))
    ax_b.add_patch(Rectangle((ix0, iy0), iw, ih, fill=False, edgecolor=TEXT_PRIMARY, linewidth=1.2))
    # the same block in the domain view
    ax_f.add_patch(Rectangle((ix0, iy0), iw, ih, fill=False, edgecolor=TEXT_PRIMARY, linewidth=1.5))
    colorbar(fig, im_b, ax_b, var)
    style(ax_b, f"block {block} with its halo\ninterior outlined, {nh[0]}-cell halo veiled")
    fig.text(0.53, 0.0, "Halo cells hold copies of the neighbouring blocks (periodic at the domain edge).",
             color=TEXT_SECONDARY, fontsize=9)

    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)


def plot_pair(a, b, labels, var, k, out):
    va, vb = var_index(a, var), var_index(b, var)
    if not np.array_equal(a["gsize"], b["gsize"]):
        sys.exit("error: the two dumps have different global sizes")
    ga, gb = global_field(a, va, k), global_field(b, vb, k)
    diff = gb - ga
    vmin, vmax = np.nanmin([ga, gb]), np.nanmax([ga, gb])

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), facecolor=SURFACE)
    for ax, d, g, label in [(axes[0], a, ga, labels[0]), (axes[1], b, gb, labels[1])]:
        im = ax.imshow(g, origin="lower", extent=domain_extent(d), cmap=SEQUENTIAL,
                       interpolation="nearest", vmin=vmin, vmax=vmax)
        draw_blocks(ax, d, k)
        colorbar(fig, im, ax, var)
        prefix = "" if label == f"step {d['step']}" else f"{label}: "
        style(ax, f"{prefix}{var} at t = {d['time']:.4g} (step {d['step']})\n{mesh_summary(d)}")

    m = np.nanmax(np.abs(diff))
    if m > 0:
        im = axes[2].imshow(diff, origin="lower", extent=domain_extent(a), cmap=DIVERGING,
                            norm=TwoSlopeNorm(0.0, -m, m), interpolation="nearest")
        colorbar(fig, im, axes[2], f"{var} difference")
        title = f"{labels[1]} - {labels[0]}, max |diff| = {m:.3e}"
    else:
        im = axes[2].imshow(np.zeros_like(diff), origin="lower", extent=domain_extent(a),
                            cmap=DIVERGING, vmin=-1, vmax=1, interpolation="nearest")
        title = f"{labels[1]} - {labels[0]}: identical (difference exactly 0)"
    style(axes[2], title)

    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=SURFACE)

    # text version: max difference of every variable
    print(f"{'variable':<10} {'max |' + labels[1] + ' - ' + labels[0] + '|':>30} {'max |' + labels[0] + '|':>16}")
    for name in a["names"]:
        if name in b["names"]:
            fa = global_field(a, a["names"].index(name), k)
            fb = global_field(b, b["names"].index(name), k)
            print(f"{name:<10} {np.nanmax(np.abs(fb - fa)):>30.3e} {np.nanmax(np.abs(fa)):>16.3e}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dumps", nargs="+", help="one or two um_dump_*.bin files")
    p.add_argument("--var", default="dens", help="variable to plot (default dens)")
    p.add_argument("--block", type=int, help="block shown with its halo (1-based; default: block with the minimum)")
    p.add_argument("--k", type=int, help="global k index of the slice in 3D (default: middle)")
    p.add_argument("--labels", nargs=2, default=None, help="labels of the two dumps")
    p.add_argument("-o", "--output", default="mesh.png", help="image file (default mesh.png)")
    args = p.parse_args()
    if len(args.dumps) > 2:
        p.error("give one or two dumps")

    dumps = [read_dump(f) for f in args.dumps]
    d0 = dumps[0]
    k = args.k if args.k is not None else (d0["gsize"][2] + 1) // 2
    print(f"{d0['path']}: {d0['ndim']}D, {mesh_summary(d0)}, t = {d0['time']:.6g}, "
          f"variables: {', '.join(d0['names'])}")

    if len(dumps) == 1:
        plot_single(d0, args.var, args.block, k, args.output)
    else:
        labels = args.labels or [f"step {d['step']}" for d in dumps]
        if labels[0] == labels[1]:
            labels = args.dumps
        plot_pair(dumps[0], dumps[1], labels, args.var, k, args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
