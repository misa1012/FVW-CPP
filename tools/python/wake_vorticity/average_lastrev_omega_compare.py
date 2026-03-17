#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
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
    p = argparse.ArgumentParser(description="Average last-revolution omega magnitude by two methods")
    p.add_argument("case_dir", help="Case directory containing wake.h5")
    p.add_argument("--config", type=str, default=None, help="Optional config json path for postprocess_grid")
    p.add_argument("--diameter", type=float, default=0.894)
    p.add_argument("--hub-z", type=float, default=0.8)
    p.add_argument("--ppd", type=float, default=100.0)
    p.add_argument("--xlim-d", nargs=2, type=float, default=[-1.0, 5.0])
    p.add_argument("--ylim-d", nargs=2, type=float, default=[-1.0, 1.0])
    p.add_argument("--zlim-d", nargs=2, type=float, default=[-1.0, 1.0])
    p.add_argument("--x-over-d-yz", type=float, default=1.0)
    p.add_argument("--steps-per-rev", type=int, default=80)
    p.add_argument("--time-stride", type=int, default=1, help="Use every Nth timestep within the selected last-revolution window")
    p.add_argument("--method", choices=["both", "circulation", "velocity"], default="both")
    p.add_argument("--output-dir", type=str, default=None, help="Optional output directory")
    p.add_argument("--timesteps", nargs="+", type=int, default=None, help="Explicit timestep list to use instead of last-revolution selection")
    return p.parse_args()


def get_case_paths(case_dir: Path):
    h5 = case_dir / "wake.h5"
    if not h5.exists():
        raise FileNotFoundError(h5)
    cfgs = sorted(case_dir.glob("config*.json"))
    if cfgs:
        return h5, cfgs[0]
    parent_cfgs = sorted(case_dir.parent.glob("config*.json"))
    case_name = case_dir.name
    matching = [p for p in parent_cfgs if case_name in p.stem]
    if matching:
        return h5, matching[0]
    if parent_cfgs:
        return h5, parent_cfgs[0]
    raise FileNotFoundError(f"No config json found in {case_dir} or {case_dir.parent}")


def get_timesteps(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        ts = sorted(int(k.split("_")[-1]) for k in f["wake"].keys() if k.startswith("timestep_"))
    return ts


def choose_last_rev_timesteps(h5_path: Path, steps_per_rev: int):
    ts = get_timesteps(h5_path)
    if len(ts) < steps_per_rev:
        raise RuntimeError("Not enough timesteps for one revolution")
    return ts[-steps_per_rev:]


def build_target_grid(D: float, hub_z: float, xlim_d, ylim_d, zlim_d, ppd: float):
    nx = int(round((xlim_d[1] - xlim_d[0]) * ppd)) + 1
    ny = int(round((ylim_d[1] - ylim_d[0]) * ppd)) + 1
    nz = int(round((zlim_d[1] - zlim_d[0]) * ppd)) + 1
    origin = np.array([xlim_d[0] * D, ylim_d[0] * D, hub_z + zlim_d[0] * D], dtype=float)
    spacing = np.array([D / ppd, D / ppd, D / ppd], dtype=float)
    dims = np.array([nx, ny, nz], dtype=int)
    upper = origin + spacing * (dims - 1)
    return origin, spacing, dims, upper


def read_lines_as_segments(h5: h5py.File, timestep: int):
    blade_ids = []
    start_xyz = []
    end_xyz = []
    tangent_xyz = []
    lengths = []
    gammas = []
    line_types = []

    grp = h5[f"/wake/timestep_{timestep}"]
    blade_names = sorted(k for k in grp.keys() if k.startswith("blade_"))
    for blade_name in blade_names:
        blade_id = int(blade_name.split("_")[-1])
        blade_grp = grp[blade_name]
        nodes = np.asarray(blade_grp["nodes"][:, :3], dtype=float)
        lines = np.asarray(blade_grp["lines"][:], dtype=float)
        if lines.size == 0:
            continue
        i0 = lines[:, 0].astype(int)
        i1 = lines[:, 1].astype(int)
        gamma = lines[:, 2].astype(float)
        ltype = lines[:, 3].astype(int)
        p0 = nodes[i0]
        p1 = nodes[i1]
        d = p1 - p0
        ds = np.linalg.norm(d, axis=1)
        good = ds > 0.0
        p0, p1, gamma, ltype, ds, d = p0[good], p1[good], gamma[good], ltype[good], ds[good], d[good]
        t_hat = d / ds[:, None]
        n = len(ds)
        blade_ids.append(np.full(n, blade_id, dtype=int))
        start_xyz.append(p0)
        end_xyz.append(p1)
        tangent_xyz.append(t_hat)
        lengths.append(ds)
        gammas.append(gamma)
        line_types.append(ltype)

    return {
        "blade_id": np.concatenate(blade_ids),
        "start_xyz": np.vstack(start_xyz),
        "end_xyz": np.vstack(end_xyz),
        "tangent_xyz": np.vstack(tangent_xyz),
        "length": np.concatenate(lengths),
        "gamma": np.concatenate(gammas),
        "line_type": np.concatenate(line_types),
    }


def subdivide_segments(start_xyz, end_xyz, vec_strength, spacing, subdivide_factor=1.0):
    h = float(np.min(spacing))
    target = subdivide_factor * h
    seg_vec = end_xyz - start_xyz
    seg_len = np.linalg.norm(seg_vec, axis=1)
    nsub = np.maximum(1, np.ceil(seg_len / target).astype(int))
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


def deposit_cic(center_xyz, vec_strength, origin, spacing, dims):
    nx, ny, nz = map(int, dims)
    omega_content = np.zeros((nx, ny, nz, 3), dtype=np.float32)
    omega_flat = omega_content.reshape(-1, 3)
    xi = (center_xyz - origin[None, :]) / spacing[None, :]
    base = np.floor(xi).astype(int)
    frac = xi - base
    valid = np.all(base >= 0, axis=1) & np.all(base < (dims - 1), axis=1)
    base = base[valid]
    frac = frac[valid]
    vec_strength = vec_strength[valid]
    flat_stride_y = nz
    flat_stride_x = ny * nz
    for ox in (0, 1):
        wx = (1.0 - frac[:, 0]) if ox == 0 else frac[:, 0]
        ix = base[:, 0] + ox
        for oy in (0, 1):
            wy = (1.0 - frac[:, 1]) if oy == 0 else frac[:, 1]
            iy = base[:, 1] + oy
            for oz in (0, 1):
                wz = (1.0 - frac[:, 2]) if oz == 0 else frac[:, 2]
                iz = base[:, 2] + oz
                w = wx * wy * wz
                flat = ix * flat_stride_x + iy * flat_stride_y + iz
                omega_flat[:, 0] += np.bincount(flat, weights=vec_strength[:, 0] * w, minlength=nx * ny * nz).astype(np.float32)
                omega_flat[:, 1] += np.bincount(flat, weights=vec_strength[:, 1] * w, minlength=nx * ny * nz).astype(np.float32)
                omega_flat[:, 2] += np.bincount(flat, weights=vec_strength[:, 2] * w, minlength=nx * ny * nz).astype(np.float32)
    return omega_content


def mean_omega_mag_from_circulation(h5_path: Path, timesteps, origin, spacing, dims):
    acc = np.zeros(tuple(dims), dtype=np.float64)
    with h5py.File(h5_path, "r") as h5:
        for i, t in enumerate(timesteps, start=1):
            print(f"[circulation] timestep {t} ({i}/{len(timesteps)})", flush=True)
            seg = read_lines_as_segments(h5, t)
            vec_strength = seg["gamma"][:, None] * seg["tangent_xyz"] * seg["length"][:, None]
            centers, strengths = subdivide_segments(seg["start_xyz"], seg["end_xyz"], vec_strength, spacing, subdivide_factor=1.0)
            omega_content = deposit_cic(centers, strengths, origin, spacing, dims)
            omega = omega_content / np.prod(spacing)
            acc += np.linalg.norm(omega, axis=3)
    return (acc / len(timesteps)).astype(np.float32)


def velocity_grid_to_omega_mag(vtk_path: Path):
    grid = pv.read(vtk_path)
    grid = grid.compute_derivative(scalars="Velocity", vorticity=True)
    dims = tuple(int(v) for v in grid.dimensions)
    vort = np.asarray(grid["vorticity"]).reshape(dims + (3,), order="F")
    mag = np.linalg.norm(vort, axis=3).astype(np.float32)
    return grid, mag


def crop_mag(grid: pv.StructuredGrid, mag: np.ndarray, D: float, hub_z: float, xlim_d, ylim_d, zlim_d):
    dims = tuple(int(v) for v in grid.dimensions)
    pts = np.asarray(grid.points).reshape(dims + (3,), order="F")
    x = pts[:, 0, 0, 0] / D
    y = pts[0, :, 0, 1] / D
    z = (pts[0, 0, :, 2] - hub_z) / D
    mx = (x >= xlim_d[0] - 1e-9) & (x <= xlim_d[1] + 1e-9)
    my = (y >= ylim_d[0] - 1e-9) & (y <= ylim_d[1] + 1e-9)
    mz = (z >= zlim_d[0] - 1e-9) & (z <= zlim_d[1] + 1e-9)
    return mag[np.ix_(mx, my, mz)]


def mean_omega_mag_from_velocity(h5_path: Path, cfg_path: Path, timesteps, D: float, hub_z: float, xlim_d, ylim_d, zlim_d, ppd: float):
    res = D / ppd
    acc = None
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"
    with tempfile.TemporaryDirectory(prefix="omega_velcurl_") as tmp:
        tmp = Path(tmp)
        for i, t in enumerate(timesteps, start=1):
            print(f"[velocity-curl] timestep {t} ({i}/{len(timesteps)})", flush=True)
            out_vtk = tmp / f"vel_t{t}.vtk"
            cmd = [
                str(POSTPROCESS_GRID),
                str(h5_path),
                str(out_vtk),
                str(cfg_path),
                str(t),
                str(res),
                "0",
                str(xlim_d[0]),
                str(xlim_d[1]),
            ]
            subprocess.run(cmd, check=True, cwd=str(PROJECT_ROOT), env=env, capture_output=True, text=True)
            grid, mag = velocity_grid_to_omega_mag(out_vtk)
            mag = crop_mag(grid, mag, D, hub_z, xlim_d, ylim_d, zlim_d)
            if acc is None:
                acc = np.zeros_like(mag, dtype=np.float64)
            acc += mag
    return (acc / len(timesteps)).astype(np.float32)


def write_image_data(path: Path, scalar: np.ndarray, origin: np.ndarray, spacing: np.ndarray, name: str):
    dims = np.array(scalar.shape, dtype=int)
    grid = pv.ImageData()
    grid.origin = tuple(origin.tolist())
    grid.spacing = tuple(spacing.tolist())
    grid.dimensions = tuple(int(v) for v in dims.tolist())
    grid.point_data[name] = scalar.reshape(-1, order="F")
    grid.save(path)


def extract_slices(arr: np.ndarray, origin: np.ndarray, spacing: np.ndarray, D: float, hub_z: float, x_over_d: float):
    nx, ny, nz = arr.shape
    x = origin[0] + spacing[0] * np.arange(nx)
    y = origin[1] + spacing[1] * np.arange(ny)
    z = origin[2] + spacing[2] * np.arange(nz)
    iz = int(np.argmin(np.abs(z - hub_z)))
    ix = int(np.argmin(np.abs(x - x_over_d * D)))
    xy = arr[:, :, iz].T
    yz = arr[ix, :, :].T
    return x / D, y / D, (z - hub_z) / D, xy, yz, x[ix] / D


def plot_slices(out_dir: Path, name: str, arr: np.ndarray, origin: np.ndarray, spacing: np.ndarray, D: float, hub_z: float, x_over_d: float):
    xD, yD, zhD, xy, yz, xD_actual = extract_slices(arr, origin, spacing, D, hub_z, x_over_d)
    vmax = float(np.nanmax(arr))

    fig, ax = plt.subplots(figsize=(7.0, 3.5), dpi=220)
    im = ax.imshow(xy, origin="lower", extent=[xD.min(), xD.max(), yD.min(), yD.max()], aspect="equal", cmap="viridis", vmin=0.0, vmax=vmax)
    ax.set_xlabel("x/D")
    ax.set_ylabel("y/D")
    ax.set_title(f"{name}: mean $|\\omega|$ on XY slice")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\overline{|\omega|}$")
    fig.tight_layout()
    fig.savefig(out_dir / f"{name}_mean_omega_mag_xy.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=220)
    im = ax.imshow(yz, origin="lower", extent=[yD.min(), yD.max(), zhD.min(), zhD.max()], aspect="equal", cmap="viridis", vmin=0.0, vmax=vmax)
    ax.set_xlabel("y/D")
    ax.set_ylabel("(z-hub)/D")
    ax.set_title(f"{name}: mean $|\\omega|$ at x/D = {xD_actual:.3f}")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\overline{|\omega|}$")
    fig.tight_layout()
    fig.savefig(out_dir / f"{name}_mean_omega_mag_yz_xD_{x_over_d:.1f}.png")
    plt.close(fig)


def main():
    args = parse_args()
    case_dir = Path(args.case_dir).resolve()
    h5_path, cfg_path = get_case_paths(case_dir)
    if args.config is not None:
        cfg_path = Path(args.config).resolve()
        if not cfg_path.exists():
            raise FileNotFoundError(cfg_path)
    out_dir = Path(args.output_dir).resolve() if args.output_dir else case_dir / "post_processing" / "omega_lastrev_compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.timesteps is not None:
        timesteps = list(args.timesteps)
    else:
        timesteps = choose_last_rev_timesteps(h5_path, args.steps_per_rev)
        timesteps = timesteps[:: args.time_stride]
    D = args.diameter
    hub_z = args.hub_z
    origin, spacing, dims, upper = build_target_grid(D, hub_z, args.xlim_d, args.ylim_d, args.zlim_d, args.ppd)

    circ_mean = None
    vel_mean = None
    if args.method in {"both", "circulation"}:
        circ_mean = mean_omega_mag_from_circulation(h5_path, timesteps, origin, spacing, dims)
        write_image_data(out_dir / "mean_omega_mag_circulation.vti", circ_mean, origin, spacing, "mean_omega_mag")
        plot_slices(out_dir, "circulation_projection", circ_mean, origin, spacing, D, hub_z, args.x_over_d_yz)
    if args.method in {"both", "velocity"}:
        vel_mean = mean_omega_mag_from_velocity(h5_path, cfg_path, timesteps, D, hub_z, args.xlim_d, args.ylim_d, args.zlim_d, args.ppd)
        write_image_data(out_dir / "mean_omega_mag_velocity_curl.vti", vel_mean, origin, spacing, "mean_omega_mag")
        plot_slices(out_dir, "velocity_curl", vel_mean, origin, spacing, D, hub_z, args.x_over_d_yz)

    if circ_mean is not None and vel_mean is not None:
        diff = circ_mean - vel_mean
        write_image_data(out_dir / "mean_omega_mag_difference_circ_minus_vel.vti", diff.astype(np.float32), origin, spacing, "mean_omega_mag_diff")
        diff_rms = float(np.sqrt(np.mean(diff**2)))
    else:
        diff_rms = None

    summary = out_dir / "summary.txt"
    summary.write_text(
        "Last-revolution mean omega magnitude comparison\n"
        f"case_dir: {case_dir}\n"
        f"timesteps: {timesteps[0]} .. {timesteps[-1]} ({len(timesteps)} snapshots)\n"
        f"D: {D}\n"
        f"hub_z: {hub_z}\n"
        f"ppd: {args.ppd}\n"
        f"domain x/D: {args.xlim_d[0]} .. {args.xlim_d[1]}\n"
        f"domain y/D: {args.ylim_d[0]} .. {args.ylim_d[1]}\n"
        f"domain (z-hub)/D: {args.zlim_d[0]} .. {args.zlim_d[1]}\n"
        f"grid dims: {dims[0]} x {dims[1]} x {dims[2]}\n"
        f"spacing: {spacing[0]:.6f}, {spacing[1]:.6f}, {spacing[2]:.6f}\n"
        + (f"circulation mean max: {float(np.max(circ_mean)):.6e}\n" if circ_mean is not None else "")
        + (f"velocity-curl mean max: {float(np.max(vel_mean)):.6e}\n" if vel_mean is not None else "")
        + (f"diff rms: {diff_rms:.6e}\n" if diff_rms is not None else "")
    )

    print(f"Saved outputs in {out_dir}")


if __name__ == "__main__":
    main()
