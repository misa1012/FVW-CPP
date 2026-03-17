#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def get_axis(origin: float, spacing: float, n: int) -> np.ndarray:
    return origin + spacing * np.arange(n)


def nearest_index(axis: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(axis - value)))


def get_field_array(grid: pv.ImageData, field: str) -> np.ndarray:
    dims = tuple(int(v) for v in grid.dimensions)
    arr = np.asarray(grid.point_data[field])
    return arr.reshape(dims, order="F")


def extract_xy(grid, field, D, hub_z, xlim_d=(-1, 5), ylim_d=(-1, 1)):
    dims = tuple(int(v) for v in grid.dimensions)
    ox, oy, oz = grid.origin
    dx, dy, dz = grid.spacing
    x = get_axis(ox, dx, dims[0]) / D
    y = get_axis(oy, dy, dims[1]) / D
    z = get_axis(oz, dz, dims[2])
    iz = nearest_index(z, hub_z)
    arr = get_field_array(grid, field)[:, :, iz].T
    mx = (x >= xlim_d[0]) & (x <= xlim_d[1])
    my = (y >= ylim_d[0]) & (y <= ylim_d[1])
    return x[mx], y[my], arr[np.ix_(my, mx)]


def extract_yz(grid, field, D, hub_z, x_over_d=1.0, lim_d=1.0):
    dims = tuple(int(v) for v in grid.dimensions)
    ox, oy, oz = grid.origin
    dx, dy, dz = grid.spacing
    x = get_axis(ox, dx, dims[0])
    y = get_axis(oy, dy, dims[1]) / D
    z = (get_axis(oz, dz, dims[2]) - hub_z) / D
    ix = nearest_index(x, x_over_d * D)
    arr = get_field_array(grid, field)[ix, :, :].T
    my = (y >= -lim_d) & (y <= lim_d)
    mz = (z >= -lim_d) & (z <= lim_d)
    return y[my], z[mz], arr[np.ix_(mz, my)]


def panel_plot(datas, titles, out, extent_kind, field, clim):
    n = len(datas)
    fig, axes = plt.subplots(1, n, figsize=(4.9 * n, 4.2), dpi=220)
    if n == 1:
        axes = [axes]
    cmap = "viridis" if field == "omega_mag" else "RdBu_r"
    im = None
    for ax, data, title in zip(axes, datas, titles):
        a, b, arr = data
        im = ax.imshow(
            arr,
            origin="lower",
            extent=[a.min(), a.max(), b.min(), b.max()],
            aspect="equal",
            cmap=cmap,
            vmin=clim[0],
            vmax=clim[1],
        )
        ax.set_title(title)
        if extent_kind == "xy":
            ax.set_xlabel("x/D")
            ax.set_ylabel("y/D")
        else:
            ax.set_xlabel("y/D")
            ax.set_ylabel("(z-hub)/D")
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.12, top=0.90, wspace=0.18)
    cax = fig.add_axes([0.92, 0.18, 0.015, 0.64])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(field)
    fig.savefig(out)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Compare 3D vorticity projection resolutions")
    ap.add_argument("base_dir", help="Directory containing ppd_75, ppd_100, ppd_150")
    ap.add_argument("--diameter", type=float, default=0.894)
    ap.add_argument("--hub-z", type=float, default=0.8)
    ap.add_argument("--x-over-d", type=float, default=1.0)
    args = ap.parse_args()

    base = Path(args.base_dir).resolve()
    ppds = [75, 100, 150]
    grids = []
    titles = []
    for ppd in ppds:
        vti = base / f"ppd_{ppd}" / "omega_grid_timestep_3200.vti"
        grids.append(pv.read(vti))
        titles.append(f"{ppd} pts/D")

    out_dir = base / "resolution_compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    xy = [extract_xy(g, "omega_z", args.diameter, args.hub_z) for g in grids]
    yz = [extract_yz(g, "omega_z", args.diameter, args.hub_z, x_over_d=args.x_over_d, lim_d=1.0) for g in grids]
    yzmag = [extract_yz(g, "omega_mag", args.diameter, args.hub_z, x_over_d=args.x_over_d, lim_d=1.0) for g in grids]

    panel_plot(xy, titles, out_dir / "omega_z_xy_resolution_compare.png", "xy", "omega_z", (-300, 300))
    panel_plot(yz, titles, out_dir / f"omega_z_yz_xD_{args.x_over_d:.1f}_resolution_compare.png", "yz", "omega_z", (-300, 300))

    maxmag = max(np.nanmax(arr) for _, _, arr in yzmag)
    panel_plot(yzmag, titles, out_dir / f"omega_mag_yz_xD_{args.x_over_d:.1f}_resolution_compare.png", "yz", "omega_mag", (0.0, maxmag))
    print(f"Saved comparison plots in {out_dir}")


if __name__ == "__main__":
    main()
