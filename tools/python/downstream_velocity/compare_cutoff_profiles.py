#!/usr/bin/env python3
"""Compare downstream velocity profiles across cutoff cases vs ALM-LES and experiment.

Generates horizontal profiles at x/D = 1,2,3,4 (default) on the hub-height plane.
FVW data are read from wake_plots/*.vtk (postprocess_grid output).
ALM-LES data are read from postProcessing/sample/<time>/horizontal_*_UMean.xy.
Experimental data are read from Markella WorkshopBT1 (TSR6_XD*_H).
"""

import argparse
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    import pyvista as pv
except Exception as e:  # pragma: no cover
    pv = None


def find_vector_name(mesh):
    # Prefer common names
    for name in ("Velocity", "velocity", "U", "u"):
        if name in mesh.array_names:
            arr = mesh[name]
            if arr.ndim == 2 and arr.shape[1] == 3:
                return name
    # Fallback: first vector
    for name in mesh.array_names:
        arr = mesh[name]
        if arr.ndim == 2 and arr.shape[1] == 3:
            return name
    raise RuntimeError("No vector array found in VTK.")


def load_fvw_vtk(vtk_path):
    if pv is None:
        raise RuntimeError("pyvista is required to read VTK files.")
    mesh = pv.read(str(vtk_path))
    vec_name = find_vector_name(mesh)
    return mesh, vec_name


def sample_horizontal_profile(mesh, vec_name, x_pos, y_min, y_max, z_hub, npts=200):
    p1 = [x_pos, y_min, z_hub]
    p2 = [x_pos, y_max, z_hub]
    line = mesh.sample_over_line(p1, p2, resolution=npts)
    pts = line.points
    U = line[vec_name]
    y = pts[:, 1]
    ux = U[:, 0]
    return y, ux


def load_alm_profiles(alm_dir, time_dir="1.872"):
    alm_dir = Path(alm_dir)
    path = alm_dir / "postProcessing" / "sample" / time_dir
    data = {}
    for f in path.glob("horizontal_*_UMean.xy"):
        label = f.name.split("_")[1]  # 1D,2D...
        arr = np.loadtxt(f)
        # columns: x y z ux uy uz
        data[f"{label}"] = {
            "y": arr[:, 1],
            "ux": arr[:, 3],
        }
    return data


def load_experiment_profiles(exp_dir):
    exp_dir = Path(exp_dir)
    data = {}
    for d in (1, 2, 3, 4):
        f = exp_dir / f"TSR6_XD{d}_H"
        if not f.exists():
            continue
        arr = np.loadtxt(f)
        # file format: z/R, (1-U/Uref), ...
        data[f"{d}D"] = {
            "y_over_r": arr[:, 0],
            "deficit": arr[:, 1],
        }
    return data


def extract_cutoff(name):
    m = re.search(r"cutoff_([0-9.]+)", name)
    return float(m.group(1)) if m else None


def main():
    parser = argparse.ArgumentParser(description="Compare cutoff study velocity profiles vs ALM/exp.")
    parser.add_argument("--study-dir", required=True, help="Cutoff study dir containing wake_plots/*.vtk")
    parser.add_argument("--alm-dir", required=True, help="ALM-LES case dir")
    parser.add_argument("--exp-dir", required=True, help="Experiment dir (WorkshopBT1)")
    parser.add_argument("--distances", nargs="+", default=["1", "2", "3", "4"],
                        help="Downstream distances in D (e.g. 1 2 3 4)")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--uinf", type=float, default=10.0)
    parser.add_argument("--r-tip", type=float, default=0.447)
    parser.add_argument("--hub-height", type=float, default=0.8)
    parser.add_argument("--y-span", type=float, default=2.0, help="Half-width in D for horizontal profile")
    parser.add_argument("--npts", type=int, default=200)
    parser.add_argument("--ylim", nargs=2, type=float, default=None, help="y-axis limits for deficit plot (e.g., 0 1)")
    args = parser.parse_args()

    study_dir = Path(args.study_dir)
    vtk_dir = study_dir / "wake_plots"
    if not vtk_dir.exists():
        raise SystemExit(f"No wake_plots found in {study_dir}")

    out_dir = Path(args.output_dir) if args.output_dir else (study_dir / "comparison_plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    alm = load_alm_profiles(args.alm_dir)
    exp = load_experiment_profiles(args.exp_dir)

    D = 2.0 * args.r_tip
    y_min = -args.y_span * D
    y_max = args.y_span * D

    # collect VTKs per cutoff
    vtk_files = sorted(vtk_dir.glob("*_cutoff_*_t*.vtk"))
    if not vtk_files:
        raise SystemExit(f"No VTK files found in {vtk_dir}")

    fvw_cases = []
    for vf in vtk_files:
        cutoff = extract_cutoff(vf.name)
        if cutoff is None:
            continue
        fvw_cases.append((cutoff, vf))

    fvw_cases.sort(key=lambda x: x[0])

    # pre-load meshes (small count) to avoid repeated IO
    meshes = []
    for cutoff, vf in fvw_cases:
        mesh, vec_name = load_fvw_vtk(vf)
        meshes.append((cutoff, mesh, vec_name))

    distances = [f"{d}D" for d in args.distances]

    for d_str in distances:
        x_pos = float(d_str.replace("D", "")) * D
        plt.figure(figsize=(6.5, 4.5))

        # FVW curves for each cutoff
        for cutoff, mesh, vec_name in meshes:
            y, ux = sample_horizontal_profile(mesh, vec_name, x_pos, y_min, y_max, args.hub_height, npts=args.npts)
            deficit = 1.0 - ux / args.uinf
            plt.plot(y / (2.0 * args.r_tip), deficit, label=rf"$\delta_c = {cutoff:.2f}$")

        # ALM (if available)
        if d_str in alm:
            y = alm[d_str]["y"]
            deficit = 1.0 - alm[d_str]["ux"] / args.uinf
            plt.plot(y / (2.0 * args.r_tip), deficit, "k--", linewidth=2, label="ALM-LES")

        # Experiment (if available)
        if d_str in exp:
            plt.plot(exp[d_str]["y_over_r"] / 2.0, exp[d_str]["deficit"], "o", color="0.6", markersize=3, label="EXP")

        plt.xlabel("y/D")
        plt.ylabel(r"$1 - U/U_{\infty}$")
        plt.title(f"Horizontal profile at x/D = {d_str}")
        plt.xlim(-1.0, 1.0)
        if args.ylim is not None:
            plt.ylim(args.ylim[0], args.ylim[1])
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=8)
        plt.tight_layout()
        out = out_dir / f"compare_horizontal_{d_str}_cutoff_study.png"
        plt.savefig(out, dpi=200)
        plt.close()
        print(f"Saved {out}")

    # 2x2 summary figure (single legend)
    if len(distances) >= 4:
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
        axes = axes.ravel()
        handles = None
        labels = None

        for idx, d_str in enumerate(distances[:4]):
            ax = axes[idx]
            x_pos = float(d_str.replace("D", "")) * D

            for cutoff, mesh, vec_name in meshes:
                y, ux = sample_horizontal_profile(mesh, vec_name, x_pos, y_min, y_max, args.hub_height, npts=args.npts)
                deficit = 1.0 - ux / args.uinf
                ax.plot(y / (2.0 * args.r_tip), deficit, label=rf"$\delta_c = {cutoff:.2f}$")

            if d_str in alm:
                y = alm[d_str]["y"]
                deficit = 1.0 - alm[d_str]["ux"] / args.uinf
                ax.plot(y / (2.0 * args.r_tip), deficit, "k--", linewidth=2, label="ALM-LES")

            if d_str in exp:
                ax.plot(exp[d_str]["y_over_r"] / 2.0, exp[d_str]["deficit"], "o", color="0.6", markersize=3, label="EXP")

            ax.set_title(f"x/D = {d_str}")
            ax.grid(True, alpha=0.3)

        for ax in axes:
            ax.set_xlim(-1.0, 1.0)
            if args.ylim is not None:
                ax.set_ylim(args.ylim[0], args.ylim[1])
            ax.set_xlabel("y/D")
        axes[0].set_ylabel(r"$1 - U/U_{\infty}$")
        axes[2].set_ylabel(r"$1 - U/U_{\infty}$")

        # single legend (use first axis handles to avoid duplicates)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, frameon=False)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        out = out_dir / "compare_horizontal_1-4D_cutoff_study.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"Saved {out}")


if __name__ == "__main__":
    main()
