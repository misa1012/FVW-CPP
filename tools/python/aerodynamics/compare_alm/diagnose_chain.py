#!/usr/bin/env python3
"""
Detailed discrepancy diagnostics for one FVW case vs one ALM case.

Outputs:
- chain_delta.png: Delta(Vrel, AoA, Cl, Cd, Fn, Ft) along span
- velocity_components_bcs.png: FVW BCS velocity components and component shares
- chain_summary.csv: spanwise values and deltas on the FVW radial grid
- chain_metrics.txt: RMSE/bias in full/mid/tip regions
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_alm_last_blade0(path: Path) -> np.ndarray:
    last = None
    with path.open("r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if int(parts[0]) == 0 and int(parts[1]) == 0:
                last = np.array([float(x) for x in parts[4:]], dtype=float)
    if last is None:
        raise RuntimeError(f"No blade0 data in {path}")
    return last


def load_fvw_spanwise(h5_path: Path, blade: int = 0, timestep: str = "last") -> dict:
    with h5py.File(h5_path, "r") as f:
        timesteps = sorted(
            int(k.split("_")[-1])
            for k in f["/liftingline"]
            if k.startswith("timestep_")
        )
        if not timesteps:
            raise RuntimeError("No timesteps found in /liftingline")
        t = timesteps[-1] if timestep == "last" else int(timestep)

        r = np.array(f["/config/geometry/r_shed"][:], dtype=float)
        r_tip = float(f["/config/simulation/r_tip"][0])

        perf = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/perf"][:], dtype=float)
        rel_bcs = np.array(
            f[f"/liftingline/timestep_{t}/blade_{blade}/relative_velocity_bcs"][:],
            dtype=float,
        )
        vmag = np.linalg.norm(rel_bcs, axis=1)

        fn = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/fn"][:], dtype=float).ravel()
        ft = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/ft"][:], dtype=float).ravel()

    return {
        "timestep": t,
        "r_norm": r / r_tip,
        "aoa": perf[:, 2],
        "cl": perf[:, 0],
        "cd": perf[:, 1],
        "fn": fn,
        "ft": ft,
        "vmag": vmag,
        "vel_bcs_x": rel_bcs[:, 0],
        "vel_bcs_y": rel_bcs[:, 1],
        "vel_bcs_z": rel_bcs[:, 2],
        "r_tip": r_tip,
    }


def load_alm_spanwise(alm_dir: Path, r_tip: float) -> dict:
    out_dir = alm_dir / "turbineOutput" / "0"

    radius = None
    with (out_dir / "radiusC").open("r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            radius = np.array([float(x) for x in line.split()[3:]], dtype=float)
            break
    if radius is None:
        raise RuntimeError("Failed to read radiusC")

    cl = read_alm_last_blade0(out_dir / "Cl")
    cd = read_alm_last_blade0(out_dir / "Cd")
    aoa = read_alm_last_blade0(out_dir / "alphaC")
    fn = read_alm_last_blade0(out_dir / "normalForce")
    ft = read_alm_last_blade0(out_dir / "tangentialForce")
    vmag = read_alm_last_blade0(out_dir / "VmagC")

    n = min(len(radius), len(cl), len(cd), len(aoa), len(fn), len(ft), len(vmag))
    radius = radius[:n]
    dr = np.gradient(radius)

    return {
        "r_norm": radius / r_tip,
        "aoa": aoa[:n],
        "cl": cl[:n],
        "cd": cd[:n],
        "fn": fn[:n] / dr,
        "ft": ft[:n] / dr,
        "vmag": vmag[:n],
    }


def plot_chain_delta(out_dir: Path, r: np.ndarray, delta: dict) -> None:
    labels = {
        "vmag": r"$\Delta |\mathbf{V}_{rel}|$ (m/s)",
        "aoa": r"$\Delta \alpha$ (deg)",
        "cl": r"$\Delta C_l$ (-)",
        "cd": r"$\Delta C_d$ (-)",
        "fn": r"$\Delta F_n$ (N/m)",
        "ft": r"$\Delta F_t$ (N/m)",
    }
    order = ["vmag", "aoa", "cl", "cd", "fn", "ft"]
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for ax, k in zip(axes.flat, order):
        ax.plot(r, delta[k], "k-o", ms=4, lw=1.2)
        ax.axhline(0.0, color="gray", lw=1.0, ls="--")
        ax.set_ylabel(labels[k])
        ax.grid(True, alpha=0.3)
    axes[2, 0].set_xlabel("r/R")
    axes[2, 1].set_xlabel("r/R")
    fig.suptitle("FVW - ALM spanwise discrepancy chain", y=0.995)
    fig.tight_layout()
    fig.savefig(out_dir / "chain_delta.png", dpi=220)
    plt.close(fig)


def plot_fvw_velocity_components(out_dir: Path, fvw: dict) -> None:
    r = fvw["r_norm"]
    vx = fvw["vel_bcs_x"]
    vy = fvw["vel_bcs_y"]
    vz = fvw["vel_bcs_z"]
    vmag = fvw["vmag"] + 1e-14

    share_x = (vx * vx) / (vmag * vmag)
    share_y = (vy * vy) / (vmag * vmag)
    share_z = (vz * vz) / (vmag * vmag)

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(r, vx, "r-o", ms=4, label=r"$V_x$ (chordwise)")
    axes[0].plot(r, vy, "b-o", ms=4, label=r"$V_y$ (normal)")
    axes[0].plot(r, vz, "g-o", ms=4, label=r"$V_z$ (spanwise)")
    axes[0].plot(r, fvw["vmag"], "k--", lw=1.2, label=r"$|\mathbf{V}_{rel}|$")
    axes[0].set_ylabel("Velocity (m/s)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="best")

    axes[1].plot(r, share_x, "r-", lw=1.6, label=r"$V_x^2/|\mathbf{V}|^2$")
    axes[1].plot(r, share_y, "b-", lw=1.6, label=r"$V_y^2/|\mathbf{V}|^2$")
    axes[1].plot(r, share_z, "g-", lw=1.6, label=r"$V_z^2/|\mathbf{V}|^2$")
    axes[1].set_ylabel("Component share (-)")
    axes[1].set_xlabel("r/R")
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="best")

    fig.suptitle("FVW relative velocity components in BCS", y=0.995)
    fig.tight_layout()
    fig.savefig(out_dir / "velocity_components_bcs.png", dpi=220)
    plt.close(fig)


def metrics_line(delta: np.ndarray, ref: np.ndarray, mask: np.ndarray) -> tuple[float, float, float]:
    d = delta[mask]
    r = ref[mask]
    bias = float(np.mean(d))
    rmse = float(np.sqrt(np.mean(d * d)))
    rel = float(100.0 * rmse / (np.sqrt(np.mean(r * r)) + 1e-12))
    return bias, rmse, rel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fvw_h5", type=Path)
    parser.add_argument("--alm-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--blade", type=int, default=0)
    parser.add_argument("--timestep", default="last")
    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir else args.fvw_h5.parent / "alm_compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    fvw = load_fvw_spanwise(args.fvw_h5, blade=args.blade, timestep=args.timestep)
    alm = load_alm_spanwise(args.alm_dir, fvw["r_tip"])

    keys = ["vmag", "aoa", "cl", "cd", "fn", "ft"]
    r = fvw["r_norm"]
    interp_alm = {}
    delta = {}
    for k in keys:
        interp_alm[k] = np.interp(r, alm["r_norm"], alm[k])
        delta[k] = fvw[k] - interp_alm[k]

    # CSV on common FVW radial grid
    csv_path = out_dir / "chain_summary.csv"
    with csv_path.open("w") as f:
        head = [
            "r_over_R",
            "fvw_vmag", "alm_vmag_interp", "delta_vmag",
            "fvw_aoa_deg", "alm_aoa_deg_interp", "delta_aoa_deg",
            "fvw_cl", "alm_cl_interp", "delta_cl",
            "fvw_cd", "alm_cd_interp", "delta_cd",
            "fvw_fn_Npm", "alm_fn_Npm_interp", "delta_fn_Npm",
            "fvw_ft_Npm", "alm_ft_Npm_interp", "delta_ft_Npm",
            "fvw_vx_bcs", "fvw_vy_bcs", "fvw_vz_bcs",
        ]
        f.write(",".join(head) + "\n")
        for i in range(len(r)):
            row = [
                r[i],
                fvw["vmag"][i], interp_alm["vmag"][i], delta["vmag"][i],
                fvw["aoa"][i], interp_alm["aoa"][i], delta["aoa"][i],
                fvw["cl"][i], interp_alm["cl"][i], delta["cl"][i],
                fvw["cd"][i], interp_alm["cd"][i], delta["cd"][i],
                fvw["fn"][i], interp_alm["fn"][i], delta["fn"][i],
                fvw["ft"][i], interp_alm["ft"][i], delta["ft"][i],
                fvw["vel_bcs_x"][i], fvw["vel_bcs_y"][i], fvw["vel_bcs_z"][i],
            ]
            f.write(",".join(f"{v:.12e}" for v in row) + "\n")

    # Region metrics
    masks = {
        "full": np.ones_like(r, dtype=bool),
        "mid_0.2_0.8": (r >= 0.2) & (r <= 0.8),
        "tip_ge_0.85": r >= 0.85,
    }
    with (out_dir / "chain_metrics.txt").open("w") as f:
        f.write("Metric format: bias, rmse, rel_rmse_percent\n")
        for k in keys:
            f.write(f"\n[{k}]\n")
            for name, m in masks.items():
                b, rmse, rel = metrics_line(delta[k], interp_alm[k], m)
                f.write(f"{name}: bias={b:.6e}, rmse={rmse:.6e}, rel={rel:.3f}%\n")
        f.write(
            "\nNote: ALM turbineOutput in this case provides |Vrel| (VmagC), "
            "but not direct per-section velocity components in blade coordinates.\n"
        )

    plot_chain_delta(out_dir, r, delta)
    plot_fvw_velocity_components(out_dir, fvw)

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {out_dir / 'chain_metrics.txt'}")
    print(f"Wrote: {out_dir / 'chain_delta.png'}")
    print(f"Wrote: {out_dir / 'velocity_components_bcs.png'}")


if __name__ == "__main__":
    main()
