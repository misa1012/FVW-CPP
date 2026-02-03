#!/usr/bin/env python3
"""
Generate VTK slices and deficit plots for cutoff study cases.
"""

import argparse
import glob
import os
import subprocess

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def run_postprocess(executable, project_root, h5_path, output_vtk, config_path,
                    timestep, resolution, slice_mode, x_start, x_end):
    cmd = [
        executable,
        h5_path,
        output_vtk,
        config_path,
        str(timestep),
        str(resolution),
        str(slice_mode),
        str(x_start),
        str(x_end),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=project_root)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "postprocess_grid failed")


def get_last_timestep(h5_path):
    with h5py.File(h5_path, "r") as f:
        if "wake" not in f:
            return 0
        steps = []
        for k in f["wake"].keys():
            if k.startswith("timestep_"):
                try:
                    steps.append(int(k.split("_")[1]))
                except Exception:
                    continue
        return max(steps) if steps else 0


def plot_deficit(vtk_path, output_png, diameter, wind_speed, case_name):
    grid = pv.read(vtk_path)
    if "Velocity" not in grid.array_names:
        raise RuntimeError("Velocity array not found in VTK")

    points = grid.points
    x = points[:, 0]
    y = points[:, 1]
    vel = grid["Velocity"]
    ux = vel[:, 0]
    deficit = (wind_speed - ux) / wind_speed

    plt.figure(figsize=(10, 4))
    levels = np.linspace(0.0, 1.0, 100)
    plt.tricontourf(x / diameter, y / diameter, deficit, levels=levels, cmap="viridis", extend="both")
    plt.colorbar(label="Deficit (1 - U/U_inf)")
    plt.xlabel("x/D")
    plt.ylabel("y/D")
    plt.title(f"Wake Deficit ({case_name})")
    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Cutoff study visualization using postprocess_grid.")
    parser.add_argument("--study-dir", required=True, help="Directory containing cutoff_* cases")
    parser.add_argument("--project-root", required=True, help="FVW-CPP project root (for data/ paths)")
    parser.add_argument("--executable", default=None, help="Path to postprocess_grid (default: <root>/build/postprocess_grid)")
    parser.add_argument("--config-template", required=True, help="Fallback config.json path")
    parser.add_argument("--timestep", type=int, default=-1, help="Target timestep (default: last)")
    parser.add_argument("--resolution", type=float, default=0.02, help="Grid resolution (m)")
    parser.add_argument("--slice-mode", type=int, default=1, help="0=3D, 1=horizontal, 2=vertical")
    parser.add_argument("--x-start", type=float, default=-1.5, help="x_start factor in D")
    parser.add_argument("--x-end", type=float, default=6.5, help="x_end factor in D")
    parser.add_argument("--diameter", type=float, default=0.894, help="Rotor diameter (m)")
    parser.add_argument("--wind-speed", type=float, default=10.0, help="Freestream wind speed (m/s)")
    parser.add_argument("--output-dir", default=None, help="Output directory for VTK and PNG")
    args = parser.parse_args()

    executable = args.executable or os.path.join(args.project_root, "build", "postprocess_grid")
    if not os.path.exists(executable):
        raise FileNotFoundError(f"Executable not found: {executable}")

    case_dirs = sorted(glob.glob(os.path.join(args.study_dir, "cutoff_*")))
    if not case_dirs:
        print(f"No cutoff_* cases found in {args.study_dir}")
        return 1

    for case_dir in case_dirs:
        case_name = os.path.basename(case_dir)
        h5_path = os.path.join(case_dir, "wake.h5")
        if not os.path.exists(h5_path):
            print(f"[{case_name}] Missing wake.h5, skipping.")
            continue

        config_path = os.path.join(case_dir, "config.json")
        if not os.path.exists(config_path):
            config_path = args.config_template

        out_dir = args.output_dir or os.path.join(case_dir, "post_processing", "vtk")
        os.makedirs(out_dir, exist_ok=True)
        vtk_path = os.path.join(out_dir, f"{case_name}_slice.vtk")
        png_path = os.path.join(out_dir, f"{case_name}_deficit.png")

        timestep = args.timestep
        if timestep < 0:
            timestep = get_last_timestep(h5_path)

        print(f"[{case_name}] Generating VTK...")
        run_postprocess(executable, args.project_root, h5_path, vtk_path, config_path,
                        timestep, args.resolution, args.slice_mode, args.x_start, args.x_end)

        print(f"[{case_name}] Plotting deficit...")
        plot_deficit(vtk_path, png_path, args.diameter, args.wind_speed, case_name)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
