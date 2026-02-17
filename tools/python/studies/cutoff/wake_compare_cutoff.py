#!/usr/bin/env python3
"""Generate 2x3 grids of velocity deficit and vorticity_z for cutoff study (last timestep)."""

import os
import json
import subprocess
import time
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv
import h5py

STUDY_DIR = Path("/home/shug8104/data/convergence_study_v2")
PROJECT_ROOT = Path("/home/shug8104/data/FVW-CPP")
EXECUTABLE = PROJECT_ROOT / "build" / "postprocess_grid"
OUT_DIR = STUDY_DIR / "wake_plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Plot settings (match visualize_fvw_wake defaults)
X_LIMITS = [-1.0, 5.0]
Y_LIMITS = [-1.0, 1.0]
VORT_VMIN, VORT_VMAX = -300, 300
DEF_VMIN, DEF_VMAX = 0.0, 1.0
PPD_VORT = 160.0


def load_diameter(model_name: str):
    params_path = PROJECT_ROOT / "data" / model_name / "turbine_params.json"
    with open(params_path) as f:
        data = json.load(f)
    r_tip = data.get("rTip", 0.447)
    return 2.0 * r_tip


def load_case_params(config_path: Path):
    with open(config_path, "r") as f:
        cfg = json.load(f)

    turb = cfg.get("turbine", {})
    sim = cfg.get("simulation", {})

    model = turb.get("model", "NTNU")
    wind_speed = turb.get("windSpeed", 10.0)
    tsr = turb.get("tsr", 6.0)
    steps_per_rev = sim.get("stepsPerRevolution", 80)

    D = load_diameter(model)
    r_tip = 0.5 * D
    omega = tsr * wind_speed / r_tip if r_tip > 1e-6 else 0.0
    dt = (2.0 * np.pi / omega) / steps_per_rev if omega > 1e-9 else 0.0

    return {"model": model, "wind_speed": wind_speed, "D": D, "dt": dt}


def last_timestep(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        timesteps = []
        if "/liftingline" in f:
            for k in f["/liftingline"].keys():
                if k.startswith("timestep_"):
                    try:
                        timesteps.append(int(k.split("_")[-1]))
                    except ValueError:
                        pass
        if not timesteps:
            # fallback: check root groups
            for k in f.keys():
                if k.startswith("timestep_"):
                    try:
                        timesteps.append(int(k.split("_")[-1]))
                    except ValueError:
                        pass
        if not timesteps:
            raise RuntimeError(f"No timesteps found in {h5_path}")
        return max(timesteps)


def run_postprocess(h5_path: Path, config_path: Path, timestep: int, resolution: float, out_vtk: Path, timeout_s: int):
    x_start_fac = X_LIMITS[0] - 0.5
    x_end_fac = X_LIMITS[1] + 0.5

    cmd = [
        str(EXECUTABLE),
        str(h5_path),
        str(out_vtk),
        str(config_path),
        str(timestep),
        str(resolution),
        "1",
        str(x_start_fac),
        str(x_end_fac),
    ]

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"

    result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(PROJECT_ROOT), timeout=timeout_s)
    if result.returncode != 0:
        raise RuntimeError(f"postprocess_grid failed for {h5_path}:\n{result.stderr}\n{result.stdout}")


def choose_grid(n):
    if n <= 1:
        return 1, 1
    if n == 2:
        return 1, 2
    if n <= 4:
        return 2, 2
    return 2, 3


def find_config_for_case(case: Path, study_dir: Path) -> Path:
    """Resolve config path for a case directory across naming variants."""
    # 1) case-local config (if user stores it inside case dir)
    local = list(case.glob("config*.json"))
    if local:
        return local[0]

    # 2) exact study-level name: config_<case_name>.json
    exact = study_dir / f"config_{case.name}.json"
    if exact.exists():
        return exact

    # 3) fallback via cutoff value in case name
    m = re.search(r"cutoff_([0-9]+(?:\\.[0-9]+)?)", case.name)
    if m:
        cutoff = m.group(1)
        candidates = [
            study_dir / f"config_cutoff_{cutoff}.json",
            study_dir / f"config_NTNU_vg_uc_cutoff_{cutoff}.json",
            study_dir / f"config_NTNU_vg_cutoff_{cutoff}.json",
        ]
        for c in candidates:
            if c.exists():
                return c

    raise RuntimeError(f"Config not found for {case}")


def case_sort_key(case: Path):
    """Sort case dirs by sweep parameter value instead of lexical name."""
    name = case.name
    m = re.search(r"_spr_(\d+)$", name)
    if m:
        return (0, float(m.group(1)), name)
    m = re.search(r"_nseg_(\d+)$", name)
    if m:
        return (1, float(m.group(1)), name)
    m = re.search(r"cutoff_([0-9]+(?:\.[0-9]+)?)", name)
    if m:
        return (2, float(m.group(1)), name)
    return (9, name)


def plot_grid(title, fields, filename, vmin, vmax, cmap, D, x_limits, y_limits, out_dir):
    n = len(fields)
    if n == 0:
        print(f"[warn] no fields to plot for {filename}")
        return

    nrows, ncols = choose_grid(n)
    fig_w = 5.0 * ncols + 1.6
    fig_h = 3.6 * nrows
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        nrows, ncols,
        left=0.06, right=0.86, bottom=0.10, top=0.92,
        wspace=0.12, hspace=0.14
    )
    axes = [fig.add_subplot(gs[i // ncols, i % ncols]) for i in range(nrows * ncols)]

    levels = np.linspace(vmin, vmax, 100)
    cf = None

    for ax, item in zip(axes, fields):
        x = item["x"] / D
        y = item["y"] / D
        val = item["val"]
        cf = ax.tricontourf(x, y, val, levels=levels, cmap=cmap, extend="both")
        ax.set_title(item["label"], fontsize=9, pad=4)
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
        ax.set_aspect("equal", "box")

    # Hide unused axes if fewer than 6 cases
    for ax in axes[len(fields):]:
        ax.set_visible(False)

    for row in range(nrows):
        ax = axes[row * ncols]
        ax.set_ylabel("y/D")
    for col in range(ncols):
        ax = axes[(nrows - 1) * ncols + col]
        ax.set_xlabel("x/D")

    # Dedicated colorbar axis to avoid overlap
    cax = fig.add_axes([0.88, 0.20, 0.015, 0.60])
    ticks = np.linspace(vmin, vmax, 5)
    cbar = fig.colorbar(cf, cax=cax, ticks=ticks)
    cbar.set_label(title)

    fig.savefig(out_dir / filename, dpi=300)
    plt.close(fig)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Wake visualization for cutoff study (2x3 grids).")
    parser.add_argument("--study-dir", default=str(STUDY_DIR), help="Directory with NTNU_cutoff_* or NTNU_*_cutoff_* cases")
    parser.add_argument("--case-dir", default=None, help="Single case directory (overrides --study-dir)")
    parser.add_argument("--ppd", type=float, default=PPD_VORT, help="Points per diameter (default: 160)")
    parser.add_argument("--res", type=float, default=None, help="Absolute grid resolution (m), overrides --ppd")
    parser.add_argument("--reuse-vtk", action="store_true", help="Reuse existing VTK if present")
    parser.add_argument("--max-cases", type=int, default=None, help="Limit number of cases (debug)")
    parser.add_argument("--xlim", nargs=2, type=float, default=X_LIMITS, help="X limits in D (e.g. -1 5)")
    parser.add_argument("--ylim", nargs=2, type=float, default=Y_LIMITS, help="Y limits in D (e.g. -1 1)")
    parser.add_argument("--timeout", type=int, default=600, help="postprocess_grid timeout (s)")
    args = parser.parse_args()

    study_dir = Path(args.study_dir)
    case_dir = Path(args.case_dir) if args.case_dir else None
    if case_dir:
        study_dir = case_dir.parent
    out_dir = (case_dir.parent / "wake_plots") if case_dir else (study_dir / "wake_plots")
    out_dir.mkdir(parents=True, exist_ok=True)
    x_limits = args.xlim
    y_limits = args.ylim
    if case_dir:
        case_dirs = [case_dir]
    else:
        case_dirs = [p for p in study_dir.glob("NTNU_cutoff_*") if p.is_dir()] + \
                    [p for p in study_dir.glob("NTNU_*_cutoff_*") if p.is_dir()]
        case_dirs = sorted(case_dirs, key=case_sort_key)
        if len(case_dirs) != 6:
            print(f"Warning: expected 6 cases, found {len(case_dirs)}")
    if args.max_cases:
        case_dirs = case_dirs[:args.max_cases]

    # load D from NTNU params (assumes same model across cases)
    D = load_diameter("NTNU")
    res_vort = args.res if args.res is not None else (D / args.ppd)

    deficit_fields = []
    vort_fields = []

    for case in case_dirs:
        h5 = case / "wake.h5"
        if not h5.exists():
            print(f"[skip] missing {h5}")
            continue

        cfg = find_config_for_case(case, study_dir)

        try:
            t_last = last_timestep(h5)
        except (OSError, RuntimeError, KeyError, ValueError) as e:
            print(f"[skip] failed to read timesteps from {h5}: {e}")
            continue
        vtk_path = out_dir / f"{case.name}_t{t_last}.vtk"
        if vtk_path.exists() and args.reuse_vtk:
            print(f"[reuse] {vtk_path}")
        else:
            print(f"[run] {case.name} t={t_last} res={res_vort:.4f}")
            t0 = time.time()
            try:
                run_postprocess(h5, cfg, t_last, res_vort, vtk_path, args.timeout)
            except (RuntimeError, OSError, subprocess.TimeoutExpired) as e:
                print(f"[skip] postprocess failed for {case.name}: {e}")
                continue
            print(f"[done] {case.name} in {time.time()-t0:.1f}s")

        try:
            grid = pv.read(vtk_path)
        except Exception as e:
            print(f"[skip] failed to read vtk {vtk_path}: {e}")
            continue
        pts = grid.points
        x = pts[:, 0]
        y = pts[:, 1]

        params = load_case_params(cfg)
        if abs(params["D"] - D) > 1e-6:
            print(f"[warn] D mismatch for {case.name}: {params['D']:.3f} vs {D:.3f}")

        # velocity deficit
        vel = grid["Velocity"]
        ux = vel[:, 0]
        wind_speed = params["wind_speed"]
        deficit = (wind_speed - ux) / wind_speed
        deficit_fields.append({
            "label": f"{case.name.replace('NTNU_', '')}, (t={t_last*params['dt']:.3f}s)",
            "x": x,
            "y": y,
            "val": deficit,
        })

        # vorticity z
        grid = grid.compute_derivative(scalars="Velocity", vorticity=True)
        vort_z = grid["vorticity"][:, 2]
        vort_fields.append({
            "label": f"{case.name.replace('NTNU_', '')}, (t={t_last*params['dt']:.3f}s)",
            "x": x,
            "y": y,
            "val": vort_z,
        })

    if not deficit_fields or not vort_fields:
        raise SystemExit("No valid cases available for wake plots.")

    plot_grid("Velocity Deficit [-]", deficit_fields, "velocity_deficit_grid.png", DEF_VMIN, DEF_VMAX, "viridis", D, x_limits, y_limits, out_dir)
    plot_grid("Vorticity $\\omega_z$ [1/s]", vort_fields, "vorticity_z_grid.png", VORT_VMIN, VORT_VMAX, "RdBu_r", D, x_limits, y_limits, out_dir)

    print("Saved:", out_dir / "velocity_deficit_grid.png")
    print("Saved:", out_dir / "vorticity_z_grid.png")


if __name__ == "__main__":
    main()
