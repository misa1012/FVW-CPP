#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


PROJECT_ROOT = Path("/home/shug8104/data/FVW-CPP")
POSTPROCESS_GRID = PROJECT_ROOT / "build" / "postprocess_grid"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Average mean |omega| on the y=0 x-z plane from velocity curl")
    p.add_argument("case_dir", help="Case directory containing wake.h5")
    p.add_argument("--config", type=str, default=None, help="Optional config json path")
    p.add_argument("--diameter", type=float, default=0.894)
    p.add_argument("--hub-z", type=float, default=0.8)
    p.add_argument("--ppd", type=float, default=60.0)
    p.add_argument("--xlim-d", nargs=2, type=float, default=[-1.0, 5.0])
    p.add_argument("--zlim-d", nargs=2, type=float, default=[-1.0, 1.0])
    p.add_argument("--timesteps", nargs="+", type=int, required=True)
    p.add_argument("--output-dir", type=str, required=True)
    return p.parse_args()


def get_case_paths(case_dir: Path):
    h5 = case_dir / "wake.h5"
    if not h5.exists():
        raise FileNotFoundError(h5)
    cfgs = sorted(case_dir.glob("config*.json"))
    if cfgs:
        return h5, cfgs[0]
    parent_cfgs = sorted(case_dir.parent.glob("config*.json"))
    matching = [p for p in parent_cfgs if case_dir.name in p.stem]
    if matching:
        return h5, matching[0]
    if parent_cfgs:
        return h5, parent_cfgs[0]
    raise FileNotFoundError(f"No config json found in {case_dir} or {case_dir.parent}")


def run_postprocess_grid(h5_path: Path, cfg_path: Path, timestep: int, grid_res: float, xlim_d, out_vtk: Path) -> None:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"
    cmd = [
        str(POSTPROCESS_GRID),
        str(h5_path),
        str(out_vtk),
        str(cfg_path),
        str(timestep),
        str(grid_res),
        "2",  # vertical plane at y=0
        str(xlim_d[0]),
        str(xlim_d[1]),
    ]
    subprocess.run(cmd, check=True, cwd=str(PROJECT_ROOT), env=env, capture_output=True, text=True)


def read_plane_omega_mag(vtk_path: Path, diameter: float, hub_z: float, zlim_d):
    grid = pv.read(vtk_path)
    grid = grid.compute_derivative(scalars="Velocity", vorticity=True)
    dims = tuple(int(v) for v in grid.dimensions)
    pts = np.asarray(grid.points).reshape(dims + (3,), order="F")
    vort = np.asarray(grid["vorticity"]).reshape(dims + (3,), order="F")
    mag = np.linalg.norm(vort, axis=3)
    xD = pts[:, 0, 0, 0] / diameter
    zD = (pts[0, 0, :, 2] - hub_z) / diameter
    mz = (zD >= zlim_d[0] - 1e-9) & (zD <= zlim_d[1] + 1e-9)
    return xD, zD[mz], mag[:, 0, mz]


def write_xz_vti(path: Path, mean_mag: np.ndarray, xD: np.ndarray, zD: np.ndarray, diameter: float, hub_z: float):
    origin = (float(xD[0] * diameter), 0.0, float(hub_z + zD[0] * diameter))
    dx = float((xD[1] - xD[0]) * diameter) if len(xD) > 1 else diameter
    dz = float((zD[1] - zD[0]) * diameter) if len(zD) > 1 else diameter
    grid = pv.ImageData()
    grid.origin = origin
    grid.spacing = (dx, 1.0, dz)
    grid.dimensions = (len(xD), 1, len(zD))
    grid.point_data["mean_omega_mag"] = mean_mag.reshape(-1, order="F").astype(np.float32)
    grid.save(path)


def plot_xz(path: Path, mean_mag: np.ndarray, xD: np.ndarray, zD: np.ndarray):
    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=240)
    im = ax.imshow(
        mean_mag.T,
        origin="lower",
        extent=[float(xD.min()), float(xD.max()), float(zD.min()), float(zD.max())],
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=float(np.nanmax(mean_mag)),
    )
    ax.set_xlabel("x/D")
    ax.set_ylabel("(z-hub)/D")
    ax.set_title(r"Mean $|\omega|$ on the $y=0$ x-z plane")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\overline{|\omega|}$")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    args = parse_args()
    case_dir = Path(args.case_dir).resolve()
    h5_path, cfg_path = get_case_paths(case_dir)
    if args.config is not None:
        cfg_path = Path(args.config).resolve()
        if not cfg_path.exists():
            raise FileNotFoundError(cfg_path)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    grid_res = args.diameter / args.ppd
    acc = None
    xD_ref = None
    zD_ref = None

    with tempfile.TemporaryDirectory(prefix="omega_velcurl_xz_") as tmp:
        tmp = Path(tmp)
        for i, timestep in enumerate(args.timesteps, start=1):
            print(f"[velocity-curl xz] timestep {timestep} ({i}/{len(args.timesteps)})", flush=True)
            out_vtk = tmp / f"vel_t{timestep}.vtk"
            run_postprocess_grid(h5_path, cfg_path, timestep, grid_res, args.xlim_d, out_vtk)
            xD, zD, mag = read_plane_omega_mag(out_vtk, args.diameter, args.hub_z, args.zlim_d)
            if acc is None:
                acc = np.zeros_like(mag, dtype=np.float64)
                xD_ref = xD
                zD_ref = zD
            acc += mag

    mean_mag = (acc / len(args.timesteps)).astype(np.float32)
    write_xz_vti(out_dir / "mean_omega_mag_velocity_curl_xz_y0.vti", mean_mag, xD_ref, zD_ref, args.diameter, args.hub_z)
    plot_xz(out_dir / "velocity_curl_mean_omega_mag_xz_y0.png", mean_mag, xD_ref, zD_ref)

    (out_dir / "summary.txt").write_text(
        "Mean |omega| from velocity curl on y=0 x-z plane\n"
        f"case_dir: {case_dir}\n"
        f"timesteps: {args.timesteps}\n"
        f"ppd: {args.ppd}\n"
        f"grid_res: {grid_res:.6f}\n"
        f"x/D: {args.xlim_d[0]} .. {args.xlim_d[1]}\n"
        f"(z-hub)/D: {args.zlim_d[0]} .. {args.zlim_d[1]}\n"
        f"grid dims: {mean_mag.shape[0]} x {mean_mag.shape[1]}\n"
        f"mean max: {float(np.max(mean_mag)):.6e}\n"
    )
    print(f"Saved outputs in {out_dir}")


if __name__ == "__main__":
    main()
