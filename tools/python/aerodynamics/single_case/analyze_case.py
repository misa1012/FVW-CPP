#!/usr/bin/env python3
"""
Analyze a single FVW case and generate selected metrics/plots.

Examples:
  conda run -n post python tools/python/aerodynamics/single_case/analyze_case.py \
      results/NTNU/wake.h5 --plots power thrust cp ct

  conda run -n post python tools/python/aerodynamics/single_case/analyze_case.py \
      results/NTNU/wake.h5 --plots aoa cl cd fn ft --timestep last
"""

import argparse
import os
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def apply_pub_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "STIXGeneral"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "axes.linewidth": 0.8,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.25,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })


def list_timesteps(f):
    ts = []
    for k in f["/liftingline"].keys():
        if k.startswith("timestep_"):
            try:
                ts.append(int(k.split("_")[-1]))
            except ValueError:
                pass
    return sorted(ts)


def compute_time_history(h5_path):
    with h5py.File(h5_path, "r") as f:
        rho = float(f["/config/simulation/rho"][0])
        U = float(f["/config/simulation/wind_speed"][0])
        r_tip = float(f["/config/simulation/r_tip"][0])
        omega = float(f["/config/simulation/omega"][0])
        dt = float(f["/config/simulation/dt"][0])
        n_b = int(f["/config/simulation/n_blades"][0])

        r_center = np.array(f["/config/geometry/r_shed"][:])

        if "/config/geometry/r_trail" in f:
            r_bounds = np.array(f["/config/geometry/r_trail"][:])
            dr = np.diff(r_bounds)
        else:
            r_hub = float(f["/config/simulation/r_hub"][0])
            dr = np.full(len(r_center), (r_tip - r_hub) / len(r_center))

        ts = list_timesteps(f)
        if not ts:
            raise RuntimeError("No timesteps found in /liftingline.")

        results = []
        for t in ts:
            torque_sum = 0.0
            thrust_sum = 0.0
            for b in range(n_b):
                blade_path = f"/liftingline/timestep_{t}/blade_{b}"
                if f"{blade_path}/fn" not in f or f"{blade_path}/ft" not in f:
                    raise RuntimeError("Missing fn/ft in wake.h5. Re-run with fn/ft output enabled.")

                fn = np.array(f[f"{blade_path}/fn"][:]).ravel()
                ft = np.array(f[f"{blade_path}/ft"][:]).ravel()

                dT = fn * dr
                dQ = ft * r_center * dr

                thrust_sum += np.sum(dT)
                torque_sum += np.sum(dQ)

            power = torque_sum * omega
            time = t * dt
            results.append([t, time, power, thrust_sum])

        results = np.array(results)

    area = np.pi * r_tip**2
    cp = results[:, 2] / (0.5 * rho * area * U**3)
    ct = results[:, 3] / (0.5 * rho * area * U**2)

    return {
        "timestep": results[:, 0],
        "time": results[:, 1],
        "power": results[:, 2],
        "thrust": results[:, 3],
        "cp": cp,
        "ct": ct,
        "omega": omega,
    }


def compute_spanwise(h5_path, timestep="last", avg_start=None, avg_end=None, blade=0, avg_blades=False):
    with h5py.File(h5_path, "r") as f:
        ts = list_timesteps(f)
        if not ts:
            raise RuntimeError("No timesteps found in /liftingline.")
        if avg_start is not None or avg_end is not None:
            start = avg_start if avg_start is not None else ts[0]
            end = avg_end if avg_end is not None else ts[-1]
            sel_ts = [t for t in ts if start <= t <= end]
            if not sel_ts:
                raise RuntimeError("No timesteps in averaging window.")
        else:
            t = ts[-1] if timestep == "last" else int(timestep)
            if t not in ts:
                raise RuntimeError(f"Timestep {t} not found.")
            sel_ts = [t]

        r = np.array(f["/config/geometry/r_shed"][:])
        r_tip = float(f["/config/simulation/r_tip"][0])
        n_b = int(f["/config/simulation/n_blades"][0])

        blades = list(range(n_b)) if avg_blades else [blade]

        aoa_acc = np.zeros_like(r, dtype=float)
        cl_acc = np.zeros_like(r, dtype=float)
        cd_acc = np.zeros_like(r, dtype=float)
        fn_acc = np.zeros_like(r, dtype=float)
        ft_acc = np.zeros_like(r, dtype=float)
        count = 0

        for t in sel_ts:
            for b in blades:
                blade_path = f"/liftingline/timestep_{t}/blade_{b}"
                perf = np.array(f[f"{blade_path}/perf"][:])
                aoa = perf[:, 2]
                cl = perf[:, 0]
                cd = perf[:, 1]

                if f"{blade_path}/fn" not in f or f"{blade_path}/ft" not in f:
                    raise RuntimeError("Missing fn/ft in wake.h5. Re-run with fn/ft output enabled.")

                fn = np.array(f[f"{blade_path}/fn"][:]).ravel()
                ft = np.array(f[f"{blade_path}/ft"][:]).ravel()

                aoa_acc += aoa
                cl_acc += cl
                cd_acc += cd
                fn_acc += fn
                ft_acc += ft
                count += 1

        if count == 0:
            raise RuntimeError("No spanwise data found.")

        return {
            "r": r,
            "r_norm": r / r_tip,
            "aoa": aoa_acc / count,
            "cl": cl_acc / count,
            "cd": cd_acc / count,
            "fn": fn_acc / count,
            "ft": ft_acc / count,
        }


def save_plot(x, y, xlabel, ylabel, title, out_path):
    apply_pub_style()
    plt.figure(figsize=(3.6, 2.6))
    plt.plot(x, y, "k-", marker="o")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Analyze a single FVW case and generate plots.")
    parser.add_argument("h5_path", help="Path to wake.h5")
    parser.add_argument("--plots", nargs="+", default=["power", "thrust", "cp", "ct"],
                        choices=["power", "thrust", "cp", "ct", "aoa", "cl", "cd", "fn", "ft"])
    parser.add_argument("--timestep", default="last", help="Spanwise timestep (int or 'last')")
    parser.add_argument("--avg-start", type=int, default=None, help="Spanwise average start timestep")
    parser.add_argument("--avg-end", type=int, default=None, help="Spanwise average end timestep")
    parser.add_argument("--blade", type=int, default=0, help="Blade index")
    parser.add_argument("--avg-blades", action="store_true", help="Average spanwise over all blades")
    parser.add_argument("--out-dir", default=None, help="Output directory")
    args = parser.parse_args()

    h5_path = Path(args.h5_path)
    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        out_dir = h5_path.parent / "aerodynamic_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    time_plots = {"power", "thrust", "cp", "ct"}
    span_plots = {"aoa", "cl", "cd", "fn", "ft"}

    if any(p in time_plots for p in args.plots):
        th = compute_time_history(str(h5_path))
        if "power" in args.plots:
            save_plot(th["time"], th["power"], "Time (s)", "Power (W)", "Power Time History",
                      out_dir / "time_power.png")
        if "thrust" in args.plots:
            save_plot(th["time"], th["thrust"], "Time (s)", "Thrust (N)", "Thrust Time History",
                      out_dir / "time_thrust.png")
        if "cp" in args.plots:
            save_plot(th["time"], th["cp"], "Time (s)", "Cp (-)", "Cp Time History",
                      out_dir / "time_cp.png")
        if "ct" in args.plots:
            save_plot(th["time"], th["ct"], "Time (s)", "Ct (-)", "Ct Time History",
                      out_dir / "time_ct.png")

        # Print last-revolution averages for Cp/Ct when requested
        if "cp" in args.plots or "ct" in args.plots:
            period = 2.0 * np.pi / th["omega"] if th["omega"] != 0.0 else None
            if period is None:
                print("Warning: omega is zero, cannot compute last-revolution averages.")
            else:
                t_end = th["time"][-1]
                t_start = t_end - period
                mask = th["time"] >= t_start
                if np.any(mask):
                    if "cp" in args.plots:
                        cp_avg = float(np.mean(th["cp"][mask]))
                        print(f"Last revolution mean Cp: {cp_avg:.4f}")
                    if "ct" in args.plots:
                        ct_avg = float(np.mean(th["ct"][mask]))
                        print(f"Last revolution mean Ct: {ct_avg:.4f}")
                else:
                    print("Warning: no samples found in last-revolution window.")

    if any(p in span_plots for p in args.plots):
        sp = compute_spanwise(str(h5_path), timestep=args.timestep,
                              avg_start=args.avg_start, avg_end=args.avg_end,
                              blade=args.blade, avg_blades=args.avg_blades)
        x = sp["r_norm"]
        if "aoa" in args.plots:
            save_plot(x, sp["aoa"], "r/R", "AoA (deg)", "AoA Spanwise", out_dir / "span_aoa.png")
        if "cl" in args.plots:
            save_plot(x, sp["cl"], "r/R", "Cl (-)", "Cl Spanwise", out_dir / "span_cl.png")
        if "cd" in args.plots:
            save_plot(x, sp["cd"], "r/R", "Cd (-)", "Cd Spanwise", out_dir / "span_cd.png")
        if "fn" in args.plots:
            save_plot(x, sp["fn"], "r/R", "Fn (N/m)", "Normal Force Spanwise", out_dir / "span_fn.png")
        if "ft" in args.plots:
            save_plot(x, sp["ft"], "r/R", "Ft (N/m)", "Tangential Force Spanwise", out_dir / "span_ft.png")

    print(f"Saved outputs to: {out_dir}")


if __name__ == "__main__":
    main()
