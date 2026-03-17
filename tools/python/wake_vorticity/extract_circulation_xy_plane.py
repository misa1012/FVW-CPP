#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract z=hub x-y plane from a 3D mean |omega| circulation field")
    p.add_argument("input_vti", help="3D VTI with point-data field mean_omega_mag")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--diameter", type=float, default=0.894)
    p.add_argument("--hub-z", type=float, default=0.8)
    p.add_argument("--vmax", type=float, default=150.0)
    return p.parse_args()


def main():
    args = parse_args()
    in_path = Path(args.input_vti).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    grid = pv.read(in_path)
    dims = tuple(int(v) for v in grid.dimensions)
    pts = np.asarray(grid.points).reshape(dims + (3,), order="F")
    arr = np.asarray(grid["mean_omega_mag"]).reshape(dims, order="F")

    zD = (pts[0, 0, :, 2] - args.hub_z) / args.diameter
    iz = int(np.argmin(np.abs(zD - 0.0)))
    xD = pts[:, 0, iz, 0] / args.diameter
    yD = pts[0, :, iz, 1] / args.diameter
    xy = arr[:, :, iz]

    # save 2D VTI
    origin = (float(xD[0] * args.diameter), float(yD[0] * args.diameter), 0.0)
    dx = float((xD[1] - xD[0]) * args.diameter) if len(xD) > 1 else args.diameter
    dy = float((yD[1] - yD[0]) * args.diameter) if len(yD) > 1 else args.diameter
    out_grid = pv.ImageData()
    out_grid.origin = origin
    out_grid.spacing = (dx, dy, 1.0)
    out_grid.dimensions = (len(xD), len(yD), 1)
    out_grid.point_data["mean_omega_mag"] = xy.reshape(-1, order="F").astype(np.float32)
    out_grid.save(out_dir / "mean_omega_mag_circulation_xy_zhub.vti")

    fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=240)
    im = ax.imshow(
        xy.T,
        origin="lower",
        extent=[float(xD.min()), float(xD.max()), float(yD.min()), float(yD.max())],
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=float(args.vmax),
    )
    ax.set_xlabel("x/D")
    ax.set_ylabel("y/D")
    ax.set_title(r"Mean $|\omega|$ on the $z=z_{\mathrm{hub}}$ x-y plane")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\overline{|\omega|}$")
    fig.tight_layout()
    fig.savefig(out_dir / "circulation_mean_omega_mag_xy_zhub.png")
    plt.close(fig)

    (out_dir / "summary.txt").write_text(
        "Mean |omega| from circulation projection on z=hub x-y plane\n"
        f"input_vti: {in_path}\n"
        f"grid dims: {xy.shape[0]} x {xy.shape[1]}\n"
        f"mean max: {float(np.max(xy)):.6e}\n"
    )

    print(f"Saved outputs in {out_dir}")


if __name__ == "__main__":
    main()
