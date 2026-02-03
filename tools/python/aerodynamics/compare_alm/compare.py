#!/usr/bin/env python3
"""
Compare FVW (wake.h5) with ALM-LES time series and spanwise quantities.

Outputs:
- time_cp_ct.png
- time_power_thrust.png
- spanwise_compare.png (AoA, Cl, Cd, Fn, Ft, |Vrel|)
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_alm_series(path):
    times = []
    values = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            times.append(float(parts[1]))
            values.append(float(parts[3]))
    return np.array(times), np.array(values)


def read_alm_last_blade0(path):
    last = None
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if int(parts[0]) == 0 and int(parts[1]) == 0:
                last = parts
    if last is None:
        raise RuntimeError(f"No blade0 data in {path}")
    return np.array([float(x) for x in last[4:]])


def load_fvw_spanwise(h5_path, blade=0, timestep="last"):
    with h5py.File(h5_path, "r") as f:
        timesteps = sorted(int(k.split("_")[-1]) for k in f["/liftingline"] if k.startswith("timestep_"))
        if not timesteps:
            raise RuntimeError("No timesteps found in /liftingline.")
        t = timesteps[-1] if timestep == "last" else int(timestep)
        if t not in timesteps:
            raise RuntimeError(f"Timestep {t} not found.")

        r = np.array(f["/config/geometry/r_shed"][:])
        r_tip = float(f["/config/simulation/r_tip"][0])

        perf = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/perf"][:])
        cl = perf[:, 0]
        cd = perf[:, 1]
        aoa = perf[:, 2]

        if f"/liftingline/timestep_{t}/blade_{blade}/fn" not in f or f"/liftingline/timestep_{t}/blade_{blade}/ft" not in f:
            raise RuntimeError("Missing fn/ft in wake.h5. Re-run with fn/ft output enabled.")

        fn = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/fn"][:]).ravel()
        ft = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/ft"][:]).ravel()

        vel = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/relative_velocity_bcs"][:])
        vmag = np.linalg.norm(vel, axis=1)

    return {
        "timestep": t,
        "r_norm": r / r_tip,
        "aoa": aoa,
        "cl": cl,
        "cd": cd,
        "fn": fn,
        "ft": ft,
        "vmag": vmag,
        "r_tip": r_tip,
    }


def load_alm_spanwise(alm_dir, r_tip):
    alm_dir = Path(alm_dir) / "turbineOutput" / "0"

    radius = None
    with open(alm_dir / "radiusC", "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            radius = np.array([float(x) for x in parts[3:]])
            break
    if radius is None:
        raise RuntimeError("Failed to read radiusC")

    cl = read_alm_last_blade0(alm_dir / "Cl")
    cd = read_alm_last_blade0(alm_dir / "Cd")
    aoa = read_alm_last_blade0(alm_dir / "alphaC")
    fn = read_alm_last_blade0(alm_dir / "normalForce")
    ft = read_alm_last_blade0(alm_dir / "tangentialForce")
    vmag = read_alm_last_blade0(alm_dir / "VmagC")

    min_len = min(len(radius), len(cl), len(cd), len(aoa), len(fn), len(ft), len(vmag))
    radius = radius[:min_len]
    cl = cl[:min_len]
    cd = cd[:min_len]
    aoa = aoa[:min_len]
    fn = fn[:min_len]
    ft = ft[:min_len]
    vmag = vmag[:min_len]

    dr = np.gradient(radius)
    fn_per_m = fn / dr
    ft_per_m = ft / dr

    return {
        "r_norm": radius / r_tip,
        "aoa": aoa,
        "cl": cl,
        "cd": cd,
        "fn": fn_per_m,
        "ft": ft_per_m,
        "vmag": vmag,
    }


def compute_fvw_time_series(h5_path):
    with h5py.File(h5_path, "r") as f:
        rho = float(f["/config/simulation/rho"][0])
        U = float(f["/config/simulation/wind_speed"][0])
        r_tip = float(f["/config/simulation/r_tip"][0])
        r_hub = float(f["/config/simulation/r_hub"][0])
        omega = float(f["/config/simulation/omega"][0])
        dt = float(f["/config/simulation/dt"][0])
        n_b = int(f["/config/simulation/n_blades"][0])

        r = np.array(f["/config/geometry/r_shed"][:])
        if "/config/geometry/r_trail" in f:
            r_bounds = np.array(f["/config/geometry/r_trail"][:])
            dr = np.diff(r_bounds)
        else:
            dr = np.full(len(r), (r_tip - r_hub) / len(r))

        timesteps = sorted(int(k.split("_")[-1]) for k in f["/liftingline"] if k.startswith("timestep_"))
        if not timesteps:
            raise RuntimeError("No timesteps found in /liftingline.")

        power = []
        thrust = []
        times = []
        for t in timesteps:
            torque_sum = 0.0
            thrust_sum = 0.0
            for b in range(n_b):
                blade_path = f"/liftingline/timestep_{t}/blade_{b}"
                if f"{blade_path}/fn" not in f or f"{blade_path}/ft" not in f:
                    raise RuntimeError("Missing fn/ft in wake.h5. Re-run with fn/ft output enabled.")
                fn = np.array(f[f"{blade_path}/fn"][:]).ravel()
                ft = np.array(f[f"{blade_path}/ft"][:]).ravel()

                thrust_sum += np.sum(fn * dr)
                torque_sum += np.sum(ft * r * dr)

            power.append(torque_sum * omega)
            thrust.append(thrust_sum)
            times.append(t * dt)

    power = np.array(power)
    thrust = np.array(thrust)
    times = np.array(times)

    area = np.pi * r_tip**2
    cp = power / (0.5 * rho * area * U**3)
    ct = thrust / (0.5 * rho * area * U**2)

    return {
        "time": times,
        "power": power,
        "thrust": thrust,
        "cp": cp,
        "ct": ct,
        "rho": rho,
        "U": U,
        "r_tip": r_tip,
    }


def plot_time_series(fvw, alm_dir, out_dir):
    alm_dir = Path(alm_dir) / "turbineOutput" / "0"
    t_p, p = read_alm_series(alm_dir / "powerRotor")
    t_t, thrust = read_alm_series(alm_dir / "thrust")

    area = np.pi * fvw["r_tip"]**2
    q_dynamic = 0.5 * fvw["rho"] * area * fvw["U"]**2
    p_dynamic = 0.5 * fvw["rho"] * area * fvw["U"]**3

    cp_alm = p / p_dynamic
    ct_alm = thrust / q_dynamic

    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax[0].plot(fvw["time"], fvw["cp"], "k-", label="FVW")
    ax[0].plot(t_p, cp_alm, "r--", label="ALM-LES")
    ax[0].set_ylabel("Cp (-)")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    ax[1].plot(fvw["time"], fvw["ct"], "k-", label="FVW")
    ax[1].plot(t_t, ct_alm, "r--", label="ALM-LES")
    ax[1].set_ylabel("Ct (-)")
    ax[1].set_xlabel("Time (s)")
    ax[1].grid(True, alpha=0.3)

    out_path = Path(out_dir) / "time_cp_ct.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax[0].plot(fvw["time"], fvw["power"], "k-", label="FVW")
    ax[0].plot(t_p, p, "r--", label="ALM-LES")
    ax[0].set_ylabel("Power (W)")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    ax[1].plot(fvw["time"], fvw["thrust"], "k-", label="FVW")
    ax[1].plot(t_t, thrust, "r--", label="ALM-LES")
    ax[1].set_ylabel("Thrust (N)")
    ax[1].set_xlabel("Time (s)")
    ax[1].grid(True, alpha=0.3)

    out_path = Path(out_dir) / "time_power_thrust.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_spanwise(fvw, alm, out_dir):
    fig, axes = plt.subplots(3, 2, figsize=(9, 9), sharex=True)

    def plot(ax, key, ylabel):
        ax.plot(fvw["r_norm"], fvw[key], "k-", marker="o", label="FVW")
        ax.plot(alm["r_norm"], alm[key], "r--", label="ALM-LES")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    plot(axes[0, 0], "aoa", "AoA (deg)")
    plot(axes[0, 1], "cl", "Cl (-)")
    plot(axes[1, 0], "cd", "Cd (-)")
    plot(axes[1, 1], "fn", "Fn (N/m)")
    plot(axes[2, 0], "ft", "Ft (N/m)")
    plot(axes[2, 1], "vmag", "|Vrel| (m/s)")

    axes[2, 0].set_xlabel("r/R")
    axes[2, 1].set_xlabel("r/R")

    axes[0, 0].legend(loc="best")

    out_path = Path(out_dir) / "spanwise_compare.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Compare FVW wake.h5 with ALM-LES outputs.")
    parser.add_argument("fvw_h5", help="Path to FVW wake.h5")
    parser.add_argument("--alm-dir", required=True, help="ALM-LES case directory (contains turbineOutput/0)")
    parser.add_argument("--out-dir", default=None, help="Output directory (default: fvw_h5 parent)")
    parser.add_argument("--timestep", default="last", help="FVW timestep for spanwise (int or 'last')")
    parser.add_argument("--blade", type=int, default=0, help="FVW blade index for spanwise")
    args = parser.parse_args()

    fvw_h5 = Path(args.fvw_h5)
    out_dir = Path(args.out_dir) if args.out_dir else fvw_h5.parent / "alm_compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    fvw_ts = compute_fvw_time_series(str(fvw_h5))
    plot_time_series(fvw_ts, args.alm_dir, out_dir)

    fvw_sp = load_fvw_spanwise(str(fvw_h5), blade=args.blade, timestep=args.timestep)
    alm_sp = load_alm_spanwise(args.alm_dir, fvw_sp["r_tip"])
    plot_spanwise(fvw_sp, alm_sp, out_dir)

    print(f"Saved comparison plots to: {out_dir}")


if __name__ == "__main__":
    main()
