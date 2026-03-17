#!/usr/bin/env python3
"""
Plot normalized XY and YZ vorticity slices from a regular-grid VTI file.

Reference plotting style:
- XY like wake_plots/vorticity_z_grid.png
- YZ like cross_section_velocity_deficit_0to3D_corrected/.../velocity_deficit_yz_*.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot normalized vorticity slices from VTI")
    p.add_argument("vti", help="Path to omega_grid_*.vti")
    p.add_argument("--out-dir", default=None, help="Output directory, default: same as VTI")
    p.add_argument("--field", default="omega_z", choices=["omega_mag", "omega_x", "omega_y", "omega_z"])
    p.add_argument("--hub-z", type=float, default=0.8, help="Hub center z [m]")
    p.add_argument("--diameter", type=float, default=0.894, help="Rotor diameter D [m]")
    p.add_argument("--x-over-d", type=float, default=1.0, help="YZ slice location x/D")
    p.add_argument("--z-xy", type=float, default=None, help="XY slice z [m], default: hub-z")
    p.add_argument("--xy-xlim", nargs=2, type=float, default=[-1.0, 5.0], help="XY x/D limits")
    p.add_argument("--xy-ylim", nargs=2, type=float, default=[-1.0, 1.0], help="XY y/D limits")
    p.add_argument("--yz-lim", type=float, default=0.75, help="YZ symmetric limit in y/D and (z-hub)/D")
    p.add_argument("--clim", nargs=2, type=float, default=None, help="Color limits")
    p.add_argument(
        "--gaussian-radius",
        type=float,
        default=0.0,
        help="If > 0, apply Gaussian smoothing in units of grid cells before plotting",
    )
    return p.parse_args()


def get_field_array(grid: pv.ImageData, field: str, gaussian_radius: float = 0.0) -> np.ndarray:
    dims = tuple(int(v) for v in grid.dimensions)
    arr = np.asarray(grid.point_data[field])
    arr = arr.reshape(dims, order="F")
    if gaussian_radius > 0.0:
        arr = gaussian_filter(arr, sigma=gaussian_radius, mode="nearest")
    return arr


def get_axis(origin: float, spacing: float, n: int) -> np.ndarray:
    return origin + spacing * np.arange(n)


def nearest_index(axis: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(axis - value)))


def plot_xy(
    grid: pv.ImageData,
    field: str,
    hub_z: float,
    D: float,
    z_xy: float | None,
    xlim_d: tuple[float, float],
    ylim_d: tuple[float, float],
    out: Path,
    clim,
    gaussian_radius: float,
) -> None:
    dims = tuple(int(v) for v in grid.dimensions)
    ox, oy, oz = grid.origin
    dx, dy, dz = grid.spacing
    x = get_axis(ox, dx, dims[0])
    y = get_axis(oy, dy, dims[1])
    z = get_axis(oz, dz, dims[2])
    z0 = hub_z if z_xy is None else z_xy
    iz = nearest_index(z, z0)

    arr3 = get_field_array(grid, field, gaussian_radius=gaussian_radius)
    arr2 = arr3[:, :, iz].T  # (ny, nx)

    xD = x / D
    yD = y / D
    mx = (xD >= xlim_d[0]) & (xD <= xlim_d[1])
    my = (yD >= ylim_d[0]) & (yD <= ylim_d[1])

    arr2 = arr2[np.ix_(my, mx)]
    xD = xD[mx]
    yD = yD[my]

    fig, ax = plt.subplots(figsize=(7.0, 3.5), dpi=220)
    im = ax.imshow(
        arr2,
        origin="lower",
        extent=[xD.min(), xD.max(), yD.min(), yD.max()],
        aspect="equal",
        cmap="viridis" if field == "omega_mag" else "RdBu_r",
        vmin=None if clim is None else clim[0],
        vmax=None if clim is None else clim[1],
    )
    ax.set_xlabel("x/D")
    ax.set_ylabel("y/D")
    ax.set_title(f"{field} on z = {z[iz]:.3f} m")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(field)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_yz(
    grid: pv.ImageData,
    field: str,
    hub_z: float,
    D: float,
    x_over_d: float,
    yz_lim: float,
    out: Path,
    clim,
    gaussian_radius: float,
) -> None:
    dims = tuple(int(v) for v in grid.dimensions)
    ox, oy, oz = grid.origin
    dx, dy, dz = grid.spacing
    x = get_axis(ox, dx, dims[0])
    y = get_axis(oy, dy, dims[1])
    z = get_axis(oz, dz, dims[2])
    x0 = x_over_d * D
    ix = nearest_index(x, x0)

    arr3 = get_field_array(grid, field, gaussian_radius=gaussian_radius)
    arr2 = arr3[ix, :, :].T  # (nz, ny)

    yD = y / D
    zhD = (z - hub_z) / D
    my = (yD >= -yz_lim) & (yD <= yz_lim)
    mz = (zhD >= -yz_lim) & (zhD <= yz_lim)

    arr2 = arr2[np.ix_(mz, my)]
    yD = yD[my]
    zhD = zhD[mz]

    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=220)
    im = ax.imshow(
        arr2,
        origin="lower",
        extent=[yD.min(), yD.max(), zhD.min(), zhD.max()],
        aspect="equal",
        cmap="viridis" if field == "omega_mag" else "RdBu_r",
        vmin=None if clim is None else clim[0],
        vmax=None if clim is None else clim[1],
    )
    ax.set_xlabel("y/D")
    ax.set_ylabel("(z - hub)/D")
    ax.set_title(f"{field} at x/D = {x[ix] / D:.3f}")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(field)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    vti = Path(args.vti).resolve()
    if not vti.exists():
        raise FileNotFoundError(vti)
    out_dir = Path(args.out_dir).resolve() if args.out_dir else vti.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    grid = pv.read(vti)
    clim = tuple(args.clim) if args.clim is not None else None

    stem = vti.stem
    out_xy = out_dir / f"{stem}_{args.field}_xy_ref.png"
    out_yz = out_dir / f"{stem}_{args.field}_yz_xD_{args.x_over_d:.1f}_ref.png"

    plot_xy(
        grid,
        args.field,
        hub_z=args.hub_z,
        D=args.diameter,
        z_xy=args.z_xy,
        xlim_d=tuple(args.xy_xlim),
        ylim_d=tuple(args.xy_ylim),
        out=out_xy,
        clim=clim,
        gaussian_radius=args.gaussian_radius,
    )
    plot_yz(
        grid,
        args.field,
        hub_z=args.hub_z,
        D=args.diameter,
        x_over_d=args.x_over_d,
        yz_lim=args.yz_lim,
        out=out_yz,
        clim=clim,
        gaussian_radius=args.gaussian_radius,
    )

    print(f"Saved {out_xy}")
    print(f"Saved {out_yz}")


if __name__ == "__main__":
    main()
