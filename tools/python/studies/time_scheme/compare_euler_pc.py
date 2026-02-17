#!/usr/bin/env python3
"""
Compare Euler vs Predictor-Corrector (PC) for one matched case.

Outputs in --out-dir:
  - wake_velocity_deficit_euler_pc.png
  - wake_vorticityz_euler_pc.png
  - spanwise_step_common.png
  - time_cp_ct_common.png
  - compare_summary.txt
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pyvista as pv

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def last_timestep(h5_path: Path) -> int:
    with h5py.File(h5_path, "r") as f:
        steps = []
        for k in f["/liftingline"].keys():
            if k.startswith("timestep_"):
                steps.append(int(k.split("_")[-1]))
        if not steps:
            raise RuntimeError(f"No /liftingline/timestep_* in {h5_path}")
        return max(steps)


def read_case_meta(config_path: Path, h5_path: Path) -> dict:
    with open(config_path, "r") as f:
        cfg = json.load(f)
    with h5py.File(h5_path, "r") as f:
        dt = float(f["/config/simulation/dt"][0])
        uinf = float(f["/config/simulation/wind_speed"][0])
        r_tip = float(f["/config/simulation/r_tip"][0])
        omega = float(f["/config/simulation/omega"][0])
    return {
        "cfg": cfg,
        "dt": dt,
        "uinf": uinf,
        "D": 2.0 * r_tip,
        "omega": omega,
    }


def run_postprocess_grid(
    executable: Path,
    project_root: Path,
    h5_path: Path,
    config_path: Path,
    timestep: int,
    res: float,
    xlim: tuple[float, float],
    out_vtk: Path,
) -> None:
    x_start_fac = xlim[0] - 0.5
    x_end_fac = xlim[1] + 0.5
    cmd = [
        str(executable),
        str(h5_path),
        str(out_vtk),
        str(config_path),
        str(timestep),
        str(res),
        "1",
        str(x_start_fac),
        str(x_end_fac),
    ]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = (
        "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/OpenSSL/1.1/lib:"
        "/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:"
        "/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:"
        "/apps/system/easybuild/software/GCCcore/11.2.0/lib64"
    )
    res_run = subprocess.run(
        cmd, capture_output=True, text=True, env=env, cwd=str(project_root)
    )
    if res_run.returncode != 0:
        raise RuntimeError(
            f"postprocess_grid failed for {h5_path}\nstdout:\n{res_run.stdout}\nstderr:\n{res_run.stderr}"
        )


def wake_fields(vtk_path: Path, uinf: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    grid = pv.read(str(vtk_path))
    pts = grid.points
    x, y = pts[:, 0], pts[:, 1]
    vel = grid["Velocity"]
    deficit = 1.0 - (vel[:, 0] / uinf)
    g2 = grid.compute_derivative(scalars="Velocity", vorticity=True)
    vortz = g2["vorticity"][:, 2]
    return x, y, deficit, vortz


def sample_euler_on_pc(
    x_e: np.ndarray, y_e: np.ndarray, f_e: np.ndarray, x_p: np.ndarray, y_p: np.ndarray
) -> np.ndarray:
    src = pv.PolyData(np.c_[x_e, y_e, np.zeros_like(x_e)])
    src["f"] = f_e
    trg = pv.PolyData(np.c_[x_p, y_p, np.zeros_like(x_p)])
    sampled = trg.sample(src)
    return np.array(sampled["f"])


def plot_pair(
    x_e: np.ndarray,
    y_e: np.ndarray,
    f_e: np.ndarray,
    x_p: np.ndarray,
    y_p: np.ndarray,
    f_p: np.ndarray,
    title: str,
    cbar: str,
    out_path: Path,
    vmin: float,
    vmax: float,
    cmap: str,
    d: float,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(10.8, 4.2), sharex=True, sharey=True, constrained_layout=True
    )
    levels = np.linspace(vmin, vmax, 120)
    cf0 = axes[0].tricontourf(x_e / d, y_e / d, f_e, levels=levels, cmap=cmap, extend="both")
    axes[1].tricontourf(x_p / d, y_p / d, f_p, levels=levels, cmap=cmap, extend="both")
    for ax, ttl in zip(axes[:2], ["Euler", "Predictor-Corrector"]):
        ax.set_title(ttl)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal", "box")
        ax.set_xlabel("x/D")
    axes[0].set_ylabel("y/D")
    c0 = fig.colorbar(cf0, ax=axes[:], shrink=0.9, pad=0.02)
    c0.set_label(cbar)

    fig.suptitle(title)
    fig.savefig(out_path, dpi=260)
    plt.close(fig)


def plot_diff(
    x_p: np.ndarray,
    y_p: np.ndarray,
    diff: np.ndarray,
    title: str,
    cbar: str,
    out_path: Path,
    d: float,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    diff_vmax: float | None = None,
) -> None:
    if diff_vmax is None:
        diff_vmax = float(np.percentile(np.abs(diff), 99.0) + 1e-12)
    fig, ax = plt.subplots(1, 1, figsize=(6.2, 4.2), constrained_layout=True)
    levels = np.linspace(-diff_vmax, diff_vmax, 120)
    cf = ax.tricontourf(x_p / d, y_p / d, diff, levels=levels, cmap="coolwarm", extend="both")
    ax.set_title("PC - Euler")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("x/D")
    ax.set_ylabel("y/D")
    cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)
    cb.set_label(cbar)
    fig.suptitle(title)
    fig.savefig(out_path, dpi=260)
    plt.close(fig)


def read_spanwise(h5_path: Path, step: int, blade: int = 0):
    with h5py.File(h5_path, "r") as f:
        r = np.array(f["/config/geometry/r_shed"][:], dtype=float)
        r_tip = float(f["/config/simulation/r_tip"][0])
        perf = np.array(f[f"/liftingline/timestep_{step}/blade_{blade}/perf"][:], dtype=float)
        fn = np.array(f[f"/liftingline/timestep_{step}/blade_{blade}/fn"][:], dtype=float).ravel()
        ft = np.array(f[f"/liftingline/timestep_{step}/blade_{blade}/ft"][:], dtype=float).ravel()
    return {"r": r / r_tip, "aoa": perf[:, 2], "cl": perf[:, 0], "fn": fn, "ft": ft}


def plot_spanwise_compare(eu: dict, pc: dict, out_path: Path) -> None:
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    items = [("aoa", "AoA (deg)"), ("cl", "Cl (-)"), ("fn", "Fn (N/m)"), ("ft", "Ft (N/m)")]
    for a, (k, ylab) in zip(ax.flat, items):
        a.plot(eu["r"], eu[k], "k-o", ms=3, label="Euler")
        a.plot(pc["r"], pc[k], "r--o", ms=3, label="PC")
        a.set_ylabel(ylab)
        a.grid(True, alpha=0.3)
    ax[1, 0].set_xlabel("r/R")
    ax[1, 1].set_xlabel("r/R")
    ax[0, 0].legend(loc="best")
    fig.suptitle("Spanwise comparison at common timestep")
    fig.tight_layout()
    fig.savefig(out_path, dpi=260)
    plt.close(fig)


def cpct_series(h5_path: Path, stride: int = 1):
    with h5py.File(h5_path, "r") as f:
        rho = float(f["/config/simulation/rho"][0])
        u = float(f["/config/simulation/wind_speed"][0])
        r_tip = float(f["/config/simulation/r_tip"][0])
        r_hub = float(f["/config/simulation/r_hub"][0])
        omega = float(f["/config/simulation/omega"][0])
        dt = float(f["/config/simulation/dt"][0])
        nb = int(f["/config/simulation/n_blades"][0])
        r = np.array(f["/config/geometry/r_shed"][:], dtype=float)
        if "/config/geometry/r_trail" in f:
            dr = np.diff(np.array(f["/config/geometry/r_trail"][:], dtype=float))
        else:
            dr = np.full(len(r), (r_tip - r_hub) / len(r))
        ts = sorted(int(k.split("_")[-1]) for k in f["/liftingline"] if k.startswith("timestep_"))
        area = np.pi * r_tip * r_tip
        cp, ct, tt = [], [], []
        ts_use = ts[::max(1, int(stride))]
        if ts_use[-1] != ts[-1]:
            ts_use.append(ts[-1])
        for t in ts_use:
            tq, th = 0.0, 0.0
            for b in range(nb):
                fn = np.array(f[f"/liftingline/timestep_{t}/blade_{b}/fn"][:], dtype=float).ravel()
                ft = np.array(f[f"/liftingline/timestep_{t}/blade_{b}/ft"][:], dtype=float).ravel()
                th += float(np.sum(fn * dr))
                tq += float(np.sum(ft * r * dr))
            cp.append((tq * omega) / (0.5 * rho * area * u**3))
            ct.append(th / (0.5 * rho * area * u**2))
            tt.append(t * dt)
    return np.array(tt), np.array(cp), np.array(ct)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--euler-case-dir", required=True)
    ap.add_argument("--pc-case-dir", required=True)
    ap.add_argument("--euler-config", required=True)
    ap.add_argument("--pc-config", required=True)
    ap.add_argument("--project-root", default="/home/shug8104/data/FVW-CPP")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--timestep", type=int, default=None, help="Force common timestep; default min(last_euler,last_pc)")
    ap.add_argument("--ppd", type=float, default=160.0)
    ap.add_argument("--xlim", nargs=2, type=float, default=[-1.0, 5.0])
    ap.add_argument("--ylim", nargs=2, type=float, default=[-1.0, 1.0])
    ap.add_argument("--def-diff-vmax", type=float, default=None, help="Color range max for deficit difference plot.")
    ap.add_argument("--vort-diff-vmax", type=float, default=None, help="Color range max for vorticity difference plot.")
    ap.add_argument("--with-time-series", action="store_true", help="Also compute and plot Cp/Ct time series (slower).")
    args = ap.parse_args()

    euler_case = Path(args.euler_case_dir)
    pc_case = Path(args.pc_case_dir)
    euler_h5 = euler_case / "wake.h5"
    pc_h5 = pc_case / "wake.h5"
    euler_cfg = Path(args.euler_config)
    pc_cfg = Path(args.pc_config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    project_root = Path(args.project_root)
    exe = project_root / "build" / "postprocess_grid"
    if not exe.exists():
        raise RuntimeError(f"Missing executable: {exe}")

    le = last_timestep(euler_h5)
    lp = last_timestep(pc_h5)
    t_common = args.timestep if args.timestep is not None else min(le, lp)

    me = read_case_meta(euler_cfg, euler_h5)
    mp = read_case_meta(pc_cfg, pc_h5)
    d = me["D"]
    res = d / args.ppd
    t_phys = t_common * me["dt"]

    vtk_e = out_dir / f"tmp_euler_t{t_common}.vtk"
    vtk_p = out_dir / f"tmp_pc_t{t_common}.vtk"
    run_postprocess_grid(exe, project_root, euler_h5, euler_cfg, t_common, res, tuple(args.xlim), vtk_e)
    run_postprocess_grid(exe, project_root, pc_h5, pc_cfg, t_common, res, tuple(args.xlim), vtk_p)

    xe, ye, defe, vortze = wake_fields(vtk_e, me["uinf"])
    xp, yp, defp, vortzp = wake_fields(vtk_p, mp["uinf"])

    plot_pair(
        xe, ye, defe, xp, yp, defp,
        title=f"Wake Velocity Deficit at step {t_common} (t={t_phys:.3f}s)",
        cbar="Velocity Deficit [-]",
        out_path=out_dir / "wake_velocity_deficit_euler_pc_1x2.png",
        vmin=0.0, vmax=1.0, cmap="viridis", d=d,
        xlim=tuple(args.xlim), ylim=tuple(args.ylim),
    )
    plot_pair(
        xe, ye, vortze, xp, yp, vortzp,
        title=f"Wake Vertical Vorticity at step {t_common} (t={t_phys:.3f}s)",
        cbar=r"Vorticity $\omega_z$ [1/s]",
        out_path=out_dir / "wake_vorticityz_euler_pc_1x2.png",
        vmin=-300.0, vmax=300.0, cmap="RdBu_r", d=d,
        xlim=tuple(args.xlim), ylim=tuple(args.ylim),
    )
    defe_diff = defp - sample_euler_on_pc(xe, ye, defe, xp, yp)
    vort_diff = vortzp - sample_euler_on_pc(xe, ye, vortze, xp, yp)
    plot_diff(
        xp, yp, defe_diff,
        title=f"Velocity Deficit Difference at step {t_common} (t={t_phys:.3f}s)",
        cbar=r"$\Delta$ Deficit [-]",
        out_path=out_dir / "wake_velocity_deficit_diff_pc_minus_euler.png",
        d=d, xlim=tuple(args.xlim), ylim=tuple(args.ylim),
        diff_vmax=args.def_diff_vmax,
    )
    plot_diff(
        xp, yp, vort_diff,
        title=f"Vorticity Difference at step {t_common} (t={t_phys:.3f}s)",
        cbar=r"$\Delta \omega_z$ [1/s]",
        out_path=out_dir / "wake_vorticityz_diff_pc_minus_euler.png",
        d=d, xlim=tuple(args.xlim), ylim=tuple(args.ylim),
        diff_vmax=args.vort_diff_vmax,
    )

    span_e = read_spanwise(euler_h5, t_common)
    span_p = read_spanwise(pc_h5, t_common)
    plot_spanwise_compare(span_e, span_p, out_dir / "spanwise_step_common.png")

    with (out_dir / "compare_summary.txt").open("w") as f:
        f.write(f"Euler last step: {le}, time={le * me['dt']:.6f}s\n")
        f.write(f"PC    last step: {lp}, time={lp * mp['dt']:.6f}s\n")
        f.write(f"Common plotted step: {t_common}, time={t_phys:.6f}s\n")

    if args.with_time_series:
        te, cpe, cte = cpct_series(euler_h5, stride=10)
        tp, cpp, ctp = cpct_series(pc_h5, stride=1)
        t_end = min(te[-1], tp[-1])
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(te[te <= t_end], cpe[te <= t_end], "k-", label="Euler")
        ax[0].plot(tp[tp <= t_end], cpp[tp <= t_end], "r--o", ms=2, label="PC")
        ax[0].set_ylabel("Cp (-)")
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(loc="best")
        ax[1].plot(te[te <= t_end], cte[te <= t_end], "k-", label="Euler")
        ax[1].plot(tp[tp <= t_end], ctp[tp <= t_end], "r--o", ms=2, label="PC")
        ax[1].set_ylabel("Ct (-)")
        ax[1].set_xlabel("Time (s)")
        ax[1].grid(True, alpha=0.3)
        fig.suptitle("Time-series comparison on common time range")
        fig.tight_layout()
        fig.savefig(out_dir / "time_cp_ct_common.png", dpi=260)
        plt.close(fig)

        trev = 2.0 * np.pi / me["omega"]
        s = t_end - trev
        m_e = (te >= s) & (te <= t_end)
        m_p = (tp >= s) & (tp <= t_end)
        with (out_dir / "compare_summary.txt").open("a") as f:
            f.write(f"Common last-rev window: [{s:.6f}, {t_end:.6f}] s\n")
            f.write(f"Cp mean (Euler): {float(np.mean(cpe[m_e])):.6f}\n")
            f.write(f"Cp mean (PC):    {float(np.mean(cpp[m_p])):.6f}\n")
            f.write(f"Ct mean (Euler): {float(np.mean(cte[m_e])):.6f}\n")
            f.write(f"Ct mean (PC):    {float(np.mean(ctp[m_p])):.6f}\n")

    print(f"Saved plots and summary to: {out_dir}")


if __name__ == "__main__":
    main()
