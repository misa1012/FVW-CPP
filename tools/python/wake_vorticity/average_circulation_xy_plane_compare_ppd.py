#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare 2D circulation-based mean |omega| on z=hub x-y plane for multiple ppd")
    p.add_argument("case_dir", help="Case directory containing wake.h5")
    p.add_argument("--diameter", type=float, default=0.894)
    p.add_argument("--hub-z", type=float, default=0.8)
    p.add_argument("--xlim-d", nargs=2, type=float, default=[-1.0, 5.0])
    p.add_argument("--ylim-d", nargs=2, type=float, default=[-1.0, 1.0])
    p.add_argument("--timesteps", nargs="+", type=int, required=True)
    p.add_argument("--ppd-list", nargs="+", type=int, default=[50, 75, 100])
    p.add_argument("--vmax", type=float, default=150.0)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def read_lines_as_segments(h5: h5py.File, timestep: int):
    start_xyz = []
    end_xyz = []
    tangent_xyz = []
    lengths = []
    gammas = []

    grp = h5[f"/wake/timestep_{timestep}"]
    blade_names = sorted(k for k in grp.keys() if k.startswith("blade_"))
    for blade_name in blade_names:
        blade_grp = grp[blade_name]
        nodes = np.asarray(blade_grp["nodes"][:, :3], dtype=float)
        lines = np.asarray(blade_grp["lines"][:], dtype=float)
        if lines.size == 0:
            continue
        i0 = lines[:, 0].astype(int)
        i1 = lines[:, 1].astype(int)
        gamma = lines[:, 2].astype(float)
        p0 = nodes[i0]
        p1 = nodes[i1]
        d = p1 - p0
        ds = np.linalg.norm(d, axis=1)
        good = ds > 0.0
        p0, p1, gamma, ds, d = p0[good], p1[good], gamma[good], ds[good], d[good]
        t_hat = d / ds[:, None]
        start_xyz.append(p0)
        end_xyz.append(p1)
        tangent_xyz.append(t_hat)
        lengths.append(ds)
        gammas.append(gamma)

    return {
        "start_xyz": np.vstack(start_xyz),
        "end_xyz": np.vstack(end_xyz),
        "tangent_xyz": np.vstack(tangent_xyz),
        "length": np.concatenate(lengths),
        "gamma": np.concatenate(gammas),
    }


def subdivide_segments(start_xyz, end_xyz, vec_strength, h):
    seg_vec = end_xyz - start_xyz
    seg_len = np.linalg.norm(seg_vec, axis=1)
    nsub = np.maximum(1, np.ceil(seg_len / h).astype(int))
    centers = []
    strengths = []
    for p0, p1, strength, ns in zip(start_xyz, end_xyz, vec_strength, nsub):
        if ns == 1:
            centers.append(0.5 * (p0 + p1))
            strengths.append(strength)
        else:
            d = (p1 - p0) / ns
            sub_strength = strength / ns
            for k in range(ns):
                centers.append(p0 + (k + 0.5) * d)
                strengths.append(sub_strength)
    return np.asarray(centers, dtype=np.float32), np.asarray(strengths, dtype=np.float32)


def deposit_xy_cic(xy: np.ndarray, vec_strength: np.ndarray, origin_xy: np.ndarray, spacing_xy: np.ndarray, dims_xy: tuple[int, int]):
    nx, ny = dims_xy
    content = np.zeros((nx, ny, 3), dtype=np.float32)
    flat = content.reshape(-1, 3)
    xi = (xy - origin_xy[None, :]) / spacing_xy[None, :]
    base = np.floor(xi).astype(int)
    frac = xi - base
    valid = np.all(base >= 0, axis=1) & (base[:, 0] < nx - 1) & (base[:, 1] < ny - 1)
    base = base[valid]
    frac = frac[valid]
    vec_strength = vec_strength[valid]

    stride_x = ny
    for ox in (0, 1):
        wx = (1.0 - frac[:, 0]) if ox == 0 else frac[:, 0]
        ix = base[:, 0] + ox
        for oy in (0, 1):
            wy = (1.0 - frac[:, 1]) if oy == 0 else frac[:, 1]
            iy = base[:, 1] + oy
            w = wx * wy
            idx = ix * stride_x + iy
            flat[:, 0] += np.bincount(idx, weights=vec_strength[:, 0] * w, minlength=nx * ny).astype(np.float32)
            flat[:, 1] += np.bincount(idx, weights=vec_strength[:, 1] * w, minlength=nx * ny).astype(np.float32)
            flat[:, 2] += np.bincount(idx, weights=vec_strength[:, 2] * w, minlength=nx * ny).astype(np.float32)
    return content


def compute_mean_xy(h5_path: Path, timesteps, D: float, hub_z: float, xlim_d, ylim_d, ppd: int):
    dx = D / ppd
    dy = dx
    dz = dx
    nx = int(round((xlim_d[1] - xlim_d[0]) * ppd)) + 1
    ny = int(round((ylim_d[1] - ylim_d[0]) * ppd)) + 1
    origin_xy = np.array([xlim_d[0] * D, ylim_d[0] * D], dtype=float)
    spacing_xy = np.array([dx, dy], dtype=float)
    xD = origin_xy[0] / D + np.arange(nx) / ppd
    yD = origin_xy[1] / D + np.arange(ny) / ppd

    acc = np.zeros((nx, ny), dtype=np.float64)
    with h5py.File(h5_path, "r") as h5:
        for i, t in enumerate(timesteps, start=1):
            print(f"[circ xy ppd={ppd}] timestep {t} ({i}/{len(timesteps)})", flush=True)
            seg = read_lines_as_segments(h5, t)
            vec_strength = seg["gamma"][:, None] * seg["tangent_xyz"] * seg["length"][:, None]
            centers, strengths = subdivide_segments(seg["start_xyz"], seg["end_xyz"], vec_strength, dx)
            mask = np.abs(centers[:, 2] - hub_z) <= 0.5 * dz
            centers = centers[mask]
            strengths = strengths[mask]
            content = deposit_xy_cic(centers[:, :2], strengths, origin_xy, spacing_xy, (nx, ny))
            omega = content / (dx * dy * dz)
            acc += np.linalg.norm(omega, axis=2)
    return xD.astype(np.float32), yD.astype(np.float32), (acc / len(timesteps)).astype(np.float32), dx


def write_vti(path: Path, arr: np.ndarray, xD: np.ndarray, yD: np.ndarray, D: float):
    origin = (float(xD[0] * D), float(yD[0] * D), 0.0)
    dx = float((xD[1] - xD[0]) * D) if len(xD) > 1 else D
    dy = float((yD[1] - yD[0]) * D) if len(yD) > 1 else D
    grid = pv.ImageData()
    grid.origin = origin
    grid.spacing = (dx, dy, 1.0)
    grid.dimensions = (len(xD), len(yD), 1)
    grid.point_data["mean_omega_mag"] = arr.reshape(-1, order="F")
    grid.save(path)


def save_plot(path: Path, arr: np.ndarray, xD: np.ndarray, yD: np.ndarray, vmax: float, title: str):
    fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=220)
    im = ax.imshow(
        arr.T,
        origin="lower",
        extent=[float(xD.min()), float(xD.max()), float(yD.min()), float(yD.max())],
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    ax.set_xlabel("x/D")
    ax.set_ylabel("y/D")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\overline{|\omega|}$")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    args = parse_args()
    case_dir = Path(args.case_dir).resolve()
    h5_path = case_dir / "wake.h5"
    if not h5_path.exists():
        raise FileNotFoundError(h5_path)
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for ppd in args.ppd_list:
        ppd_dir = out_dir / f"ppd_{ppd}"
        ppd_dir.mkdir(exist_ok=True)
        xD, yD, mean_mag, dx = compute_mean_xy(h5_path, args.timesteps, args.diameter, args.hub_z, args.xlim_d, args.ylim_d, ppd)
        write_vti(ppd_dir / f"mean_omega_mag_circulation_xy_zhub_ppd{ppd}.vti", mean_mag, xD, yD, args.diameter)
        save_plot(
            ppd_dir / f"circulation_mean_omega_mag_xy_zhub_ppd{ppd}.png",
            mean_mag, xD, yD, args.vmax,
            rf"Circulation mean $|\omega|$ on $z=z_{{\mathrm{{hub}}}}$, ppd={ppd}",
        )
        (ppd_dir / "summary.txt").write_text(
            "2D circulation mean |omega| on z=hub x-y plane\n"
            f"timesteps: {args.timesteps}\n"
            f"ppd: {ppd}\n"
            f"grid dims: {mean_mag.shape[0]} x {mean_mag.shape[1]}\n"
            f"grid_res: {dx:.6f}\n"
            f"mean max: {float(np.max(mean_mag)):.6e}\n"
        )
        results[ppd] = (xD, yD, mean_mag)

    # comparison panel
    fig, axes = plt.subplots(1, len(args.ppd_list), figsize=(4.8 * len(args.ppd_list), 4.0), dpi=220, constrained_layout=True)
    if len(args.ppd_list) == 1:
        axes = [axes]
    for ax, ppd in zip(axes, args.ppd_list):
        xD, yD, mean_mag = results[ppd]
        im = ax.imshow(
            mean_mag.T,
            origin="lower",
            extent=[float(xD.min()), float(xD.max()), float(yD.min()), float(yD.max())],
            aspect="equal",
            cmap="viridis",
            vmin=0.0,
            vmax=args.vmax,
        )
        ax.set_title(f"ppd = {ppd}")
        ax.set_xlabel("x/D")
        ax.set_ylabel("y/D")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$\overline{|\omega|}$")
    fig.savefig(out_dir / "compare_ppd_50_75_100.png")
    plt.close(fig)

    print(f"Saved outputs in {out_dir}")


if __name__ == "__main__":
    main()
