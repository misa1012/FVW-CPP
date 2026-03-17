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
    p = argparse.ArgumentParser(description="Average mean |omega| on a blade-fixed rotating plane from velocity curl")
    p.add_argument("case_dir", help="Case directory containing wake.h5")
    p.add_argument("--config", type=str, default=None, help="Optional config json path")
    p.add_argument("--diameter", type=float, default=0.894)
    p.add_argument("--hub-z", type=float, default=0.8)
    p.add_argument("--ppd", type=float, default=100.0)
    p.add_argument("--xlim-d", nargs=2, type=float, default=[-1.0, 5.0])
    p.add_argument("--s-span-d", type=float, default=1.5, help="Half span in plane coordinate s/D")
    p.add_argument("--timesteps", nargs="+", type=int, required=True)
    p.add_argument("--blade-id", type=int, default=0)
    p.add_argument("--vmax", type=float, default=150.0)
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


def blade_plane_angle_deg(h5: h5py.File, timestep: int, blade_id: int, hub_z: float) -> float:
    grp = h5[f"/wake/timestep_{timestep}/blade_{blade_id}"]
    nodes = np.asarray(grp["nodes"][:, :3], dtype=float)
    lines = np.asarray(grp["lines"][:], dtype=float)
    bound = lines[lines[:, 3].astype(int) == 0]
    if bound.size == 0:
        raise RuntimeError(f"No bound lines found for blade {blade_id} at timestep {timestep}")
    idx = np.unique(bound[:, :2].astype(int).ravel())
    bnodes = nodes[idx]
    radial = np.sqrt(bnodes[:, 1] ** 2 + (bnodes[:, 2] - hub_z) ** 2)
    tip = bnodes[np.argmax(radial)]
    return float(np.degrees(np.arctan2(tip[1], tip[2] - hub_z)))


def run_postprocess_grid(h5_path: Path, cfg_path: Path, timestep: int, grid_res: float, xlim_d, plane_angle_deg: float, out_vtk: Path) -> None:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"
    cmd = [
        str(POSTPROCESS_GRID),
        str(h5_path),
        str(out_vtk),
        str(cfg_path),
        str(timestep),
        str(grid_res),
        "3",  # rotated vertical plane
        str(xlim_d[0]),
        str(xlim_d[1]),
        str(plane_angle_deg),
    ]
    subprocess.run(cmd, check=True, cwd=str(PROJECT_ROOT), env=env, capture_output=True, text=True)


def read_plane_omega_mag(vtk_path: Path, diameter: float, hub_z: float, alpha_deg: float, s_span_d: float, ppd: float):
    grid = pv.read(vtk_path)
    grid = grid.compute_derivative(scalars="Velocity", vorticity=True)
    dims = tuple(int(v) for v in grid.dimensions)
    pts = np.asarray(grid.points).reshape(dims + (3,), order="F")
    vort = np.asarray(grid["vorticity"]).reshape(dims + (3,), order="F")
    mag = np.linalg.norm(vort, axis=3)

    alpha = np.radians(alpha_deg)
    xD = pts[:, 0, 0, 0] / diameter
    sD = (pts[0, :, 0, 1] * np.sin(alpha) + (pts[0, :, 0, 2] - hub_z) * np.cos(alpha)) / diameter
    ms = (sD >= -s_span_d - 1e-4) & (sD <= s_span_d + 1e-4)
    s_sel = sD[ms]
    mag_sel = mag[:, ms, 0]
    s_target = -s_span_d + np.arange(int(round(2.0 * s_span_d * ppd)) + 1, dtype=np.float32) / ppd
    mag_target = np.empty((mag_sel.shape[0], s_target.size), dtype=np.float32)
    for i in range(mag_sel.shape[0]):
        mag_target[i, :] = np.interp(s_target, s_sel, mag_sel[i, :]).astype(np.float32)
    return xD, s_target, mag_target


def write_plane_vti(path: Path, mean_mag: np.ndarray, xD: np.ndarray, sD: np.ndarray, diameter: float):
    origin = (float(xD[0] * diameter), float(sD[0] * diameter), 0.0)
    dx = float((xD[1] - xD[0]) * diameter) if len(xD) > 1 else diameter
    ds = float((sD[1] - sD[0]) * diameter) if len(sD) > 1 else diameter
    grid = pv.ImageData()
    grid.origin = origin
    grid.spacing = (dx, ds, 1.0)
    grid.dimensions = (len(xD), len(sD), 1)
    grid.point_data["mean_omega_mag"] = mean_mag.reshape(-1, order="F").astype(np.float32)
    grid.save(path)


def plot_plane(path: Path, mean_mag: np.ndarray, xD: np.ndarray, sD: np.ndarray, vmax: float):
    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=240)
    im = ax.imshow(
        mean_mag.T,
        origin="lower",
        extent=[float(xD.min()), float(xD.max()), float(sD.min()), float(sD.max())],
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=float(vmax),
    )
    ax.set_xlabel("x/D")
    ax.set_ylabel("s/D")
    ax.set_title(r"Mean $|\omega|$ on a blade-fixed plane")
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

    with h5py.File(h5_path, "r") as h5:
        plane_angles = [blade_plane_angle_deg(h5, t, args.blade_id, args.hub_z) for t in args.timesteps]

    grid_res = args.diameter / args.ppd
    acc = None
    xD_ref = None
    sD_ref = None

    with tempfile.TemporaryDirectory(prefix="omega_bladefixed_") as tmp:
        tmp = Path(tmp)
        for i, (timestep, angle_deg) in enumerate(zip(args.timesteps, plane_angles), start=1):
            print(f"[blade-fixed velocity-curl] timestep {timestep} ({i}/{len(args.timesteps)}), alpha={angle_deg:.3f} deg", flush=True)
            out_vtk = tmp / f"vel_t{timestep}.vtk"
            run_postprocess_grid(h5_path, cfg_path, timestep, grid_res, args.xlim_d, angle_deg, out_vtk)
            xD, sD, mag = read_plane_omega_mag(out_vtk, args.diameter, args.hub_z, angle_deg, args.s_span_d, args.ppd)
            if acc is None:
                acc = np.zeros_like(mag, dtype=np.float64)
                xD_ref = xD
                sD_ref = sD
            acc += mag

    mean_mag = (acc / len(args.timesteps)).astype(np.float32)
    write_plane_vti(out_dir / "mean_omega_mag_velocity_curl_bladefixed_plane.vti", mean_mag, xD_ref, sD_ref, args.diameter)
    plot_plane(out_dir / "velocity_curl_mean_omega_mag_bladefixed_plane.png", mean_mag, xD_ref, sD_ref, args.vmax)
    (out_dir / "summary.txt").write_text(
        "Mean |omega| from velocity curl on a blade-fixed rotating plane\n"
        f"case_dir: {case_dir}\n"
        f"blade_id: {args.blade_id}\n"
        f"timesteps: {args.timesteps}\n"
        f"plane angles deg: {[round(a, 6) for a in plane_angles]}\n"
        f"ppd: {args.ppd}\n"
        f"grid_res: {grid_res:.6f}\n"
        f"x/D: {args.xlim_d[0]} .. {args.xlim_d[1]}\n"
        f"s/D: {-args.s_span_d} .. {args.s_span_d}\n"
        f"grid dims: {mean_mag.shape[0]} x {mean_mag.shape[1]}\n"
        f"mean max: {float(np.max(mean_mag)):.6e}\n"
    )
    print(f"Saved outputs in {out_dir}")


if __name__ == "__main__":
    main()
