#!/usr/bin/env python3
"""Postprocess cutoff study: compare Cp/Ct vs cutoff and ALM baseline.

Usage examples:
  # Default: cp/ct + spanwise AoA/Cl/Fn/Ft vs ALM
  conda run -n post python postprocess_cutoff.py

  # Only Cp/Ct convergence curves
  conda run -n post python postprocess_cutoff.py --plots cp ct

  # Spanwise only (AoA/Cl/Fn/Ft) with ALM overlay
  conda run -n post python postprocess_cutoff.py --plots aoa cl fn ft

Notes:
  - Expects case dirs under --study-dir named NTNU_cutoff_0.xx or NTNU_*_cutoff_0.xx
  - For a single case, use --case-dir /path/to/NTNU_*_cutoff_0.xx
  - Each case dir must contain wake.h5 with fn/ft datasets
  - ALM baseline is read from --alm-dir (turbineOutput/0)
"""

from pathlib import Path
import re
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY_DIR = Path("/home/shug8104/data/convergence_study")
ALM_DIR = Path("/home/shug8104/data/ALM-LES/03_v2_Refined_NoNacelle_nx150_nt240")


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

    aoa = read_alm_last_blade0(alm_dir / "alphaC")
    cl = read_alm_last_blade0(alm_dir / "Cl")
    fn = read_alm_last_blade0(alm_dir / "normalForce")
    ft = read_alm_last_blade0(alm_dir / "tangentialForce")
    vmag = read_alm_last_blade0(alm_dir / "VmagC")

    min_len = min(len(radius), len(aoa), len(cl), len(fn), len(ft), len(vmag))
    radius = radius[:min_len]
    aoa = aoa[:min_len]
    cl = cl[:min_len]
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
        "fn": fn_per_m,
        "ft": ft_per_m,
        "vmag": vmag,
    }


def compute_time_series(h5_path: Path):
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
                    raise RuntimeError(f"Missing fn/ft in {h5_path}")
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
        "omega": omega,
    }


def compute_spanwise(h5_path: Path, blade: int = 0, timestep: str = "last"):
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
        aoa = perf[:, 2]

        if f"/liftingline/timestep_{t}/blade_{blade}/fn" not in f or f"/liftingline/timestep_{t}/blade_{blade}/ft" not in f:
            raise RuntimeError("Missing fn/ft in wake.h5")
        fn = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/fn"][:]).ravel()
        ft = np.array(f[f"/liftingline/timestep_{t}/blade_{blade}/ft"][:]).ravel()

    return {
        "timestep": t,
        "r_norm": r / r_tip,
        "aoa": aoa,
        "cl": cl,
        "fn": fn,
        "ft": ft,
    }


def last_rev_mean(time, values, omega):
    period = 2.0 * np.pi / omega
    t_end = time[-1]
    t_start = t_end - period
    mask = time >= t_start
    return float(np.mean(values[mask])) if np.any(mask) else float(np.mean(values))


def parse_study_param(case_name: str):
    """Parse sweep parameter from case name.

    Supported patterns:
      - ..._nseg_<int>
      - ..._spr_<int>
      - ..._cutoff_<float>
    """
    m = re.search(r"_nseg_(\d+)$", case_name)
    if m:
        v = int(m.group(1))
        return {
            "key": "nseg",
            "value": float(v),
            "x_label": r"$N_{\mathrm{seg}}$",
            "legend_title": r"$N_{\mathrm{seg}}$",
            "label": rf"$N_{{\mathrm{{seg}}}}={v}$",
        }

    m = re.search(r"_spr_(\d+)$", case_name)
    if m:
        v = int(m.group(1))
        return {
            "key": "spr",
            "value": float(v),
            "x_label": r"$N_{\mathrm{spr}}$",
            "legend_title": r"$N_{\mathrm{spr}}$",
            "label": rf"$N_{{\mathrm{{spr}}}}={v}$",
        }

    m = re.search(r"cutoff_([0-9]+(?:\.[0-9]+)?)", case_name)
    if m:
        v = float(m.group(1))
        return {
            "key": "cutoff",
            "value": v,
            "x_label": r"$\delta_c$",
            "legend_title": r"$\delta_c$",
            "label": rf"$\delta_c={v:.2f}$",
        }

    return None


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Postprocess cutoff convergence study.")
    parser.add_argument("--study-dir", default=str(STUDY_DIR), help="Directory containing NTNU_cutoff_* or NTNU_*_cutoff_* cases")
    parser.add_argument("--case-dir", default=None, help="Single case directory (overrides --study-dir)")
    parser.add_argument("--alm-dir", default=str(ALM_DIR), help="ALM-LES case directory")
    parser.add_argument("--out-dir", default=None, help="Output directory (default: <study>/postprocess)")
    parser.add_argument("--plots", nargs="+", default=["cp", "ct", "aoa", "cl", "fn", "ft"],
                        help="Plots to generate: cp ct aoa cl fn ft")
    args = parser.parse_args()

    study_dir = Path(args.study_dir)
    case_dir = Path(args.case_dir) if args.case_dir else None
    out_dir = Path(args.out_dir) if args.out_dir else (case_dir / "postprocess" if case_dir else study_dir / "postprocess")
    out_dir.mkdir(parents=True, exist_ok=True)

    if case_dir:
        case_dirs = [case_dir]
    else:
        case_dirs = sorted([p for p in study_dir.glob("NTNU_cutoff_*") if p.is_dir()] +
                           [p for p in study_dir.glob("NTNU_*_cutoff_*") if p.is_dir()])
        if not case_dirs:
            raise SystemExit("No NTNU_cutoff_* or NTNU_*_cutoff_* case dirs found.")

    summary = []
    spanwise = []
    ref_params = None
    study_param_meta = None

    for case in case_dirs:
        h5 = case / "wake.h5"
        if not h5.exists():
            print(f"[skip] missing {h5}")
            continue
        try:
            ts = compute_time_series(h5)
        except (OSError, RuntimeError, KeyError, ValueError) as e:
            print(f"[skip] failed to read {h5}: {e}")
            continue
        p = parse_study_param(case.name)
        if p is None:
            print(f"[skip] cannot parse sweep parameter from case name: {case.name}")
            continue
        if study_param_meta is None:
            study_param_meta = p
        elif p["key"] != study_param_meta["key"]:
            print(f"[skip] mixed sweep types not supported in one run: {case.name}")
            continue
        cp_mean = last_rev_mean(ts["time"], ts["cp"], ts["omega"])
        ct_mean = last_rev_mean(ts["time"], ts["ct"], ts["omega"])
        power_mean = last_rev_mean(ts["time"], ts["power"], ts["omega"])
        thrust_mean = last_rev_mean(ts["time"], ts["thrust"], ts["omega"])

        if ref_params is None:
            ref_params = ts

        summary.append([p["value"], cp_mean, ct_mean, power_mean, thrust_mean])

        try:
            sp = compute_spanwise(h5)
        except (OSError, RuntimeError, KeyError, ValueError) as e:
            print(f"[skip] failed spanwise read {h5}: {e}")
            continue
        tip_idx = int(np.argmax(sp["r_norm"]))
        spanwise.append({
            "param": p["value"],
            "param_label": p["label"],
            "r_norm": sp["r_norm"],
            "aoa": sp["aoa"],
            "cl": sp["cl"],
            "fn": sp["fn"],
            "ft": sp["ft"],
            "fn_tip": float(sp["fn"][tip_idx]),
            "ft_tip": float(sp["ft"][tip_idx]),
        })

    summary = sorted(summary, key=lambda x: x[0])
    spanwise = sorted(spanwise, key=lambda x: x["param"])

    if not summary:
        raise SystemExit("No valid cases found (all missing/corrupt/incomplete).")
    if study_param_meta is None:
        raise SystemExit("No valid cases found with parseable sweep parameter.")

    # ALM baseline
    alm_dir = ALM_DIR / "turbineOutput" / "0"
    t_p, p = read_alm_series(alm_dir / "powerRotor")
    t_t, thrust = read_alm_series(alm_dir / "thrust")

    area = np.pi * ref_params["r_tip"]**2
    q_dynamic = 0.5 * ref_params["rho"] * area * ref_params["U"]**2
    p_dynamic = 0.5 * ref_params["rho"] * area * ref_params["U"]**3

    cp_alm = p / p_dynamic
    ct_alm = thrust / q_dynamic
    cp_alm_mean = last_rev_mean(t_p, cp_alm, ref_params["omega"])
    ct_alm_mean = last_rev_mean(t_t, ct_alm, ref_params["omega"])

    # Save CSV
    csv_path = out_dir / "cutoff_summary.csv"
    with open(csv_path, "w") as f:
        f.write(f"{study_param_meta['key']},cp_mean,ct_mean,power_mean,thrust_mean\n")
        for row in summary:
            f.write(",".join(f"{v:.6f}" for v in row) + "\n")
    print("Wrote", csv_path)

    # Plot Cp/Ct vs cutoff
    cutoffs = [r[0] for r in summary]
    cp_means = [r[1] for r in summary]
    ct_means = [r[2] for r in summary]

    plt.figure(figsize=(6.5, 4.5))
    plt.plot(cutoffs, cp_means, "o-", label="FVW Cp")
    plt.axhline(cp_alm_mean, color="r", linestyle="--", label="ALM Cp")
    plt.xlabel("cutoffParam")
    plt.ylabel("Cp")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    if "cp" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        plt.plot(cutoffs, cp_means, "o-", label="FVW Cp")
        plt.axhline(cp_alm_mean, color="r", linestyle="--", label="ALM Cp")
        plt.xlabel(study_param_meta["x_label"])
        plt.ylabel("Cp")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "cp_vs_cutoff.png", dpi=200)
        plt.close()

    if "ct" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        plt.plot(cutoffs, ct_means, "o-", label="FVW Ct")
        plt.axhline(ct_alm_mean, color="r", linestyle="--", label="ALM Ct")
        plt.xlabel(study_param_meta["x_label"])
        plt.ylabel("Ct")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "ct_vs_cutoff.png", dpi=200)
        plt.close()

    # Tip Fn/Ft vs cutoff
    _ = [r["param"] for r in spanwise]

    # AoA/Cl/Fn/Ft spanwise by cutoff
    alm_sp = None
    if any(p in args.plots for p in ["aoa", "cl", "fn", "ft"]):
        alm_sp = load_alm_spanwise(args.alm_dir, ref_params["r_tip"])

    if "aoa" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        for r in spanwise:
            plt.plot(r["r_norm"], r["aoa"], ".-", label=r["param_label"])
        if alm_sp is not None:
            plt.plot(alm_sp["r_norm"], alm_sp["aoa"], "k--", label="ALM-LES")
        plt.xlabel("r/R")
        plt.ylabel("AoA (deg)")
        plt.grid(True, alpha=0.3)
        plt.legend(title=study_param_meta["legend_title"], fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / "aoa_by_cutoff.png", dpi=200)
        plt.close()

    if "cl" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        for r in spanwise:
            plt.plot(r["r_norm"], r["cl"], ".-", label=r["param_label"])
        if alm_sp is not None:
            plt.plot(alm_sp["r_norm"], alm_sp["cl"], "k--", label="ALM-LES")
        plt.xlabel("r/R")
        plt.ylabel("Cl (-)")
        plt.grid(True, alpha=0.3)
        plt.legend(title=study_param_meta["legend_title"], fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / "cl_by_cutoff.png", dpi=200)
        plt.close()

    if "fn" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        for r in spanwise:
            plt.plot(r["r_norm"], r["fn"], ".-", label=r["param_label"])
        if alm_sp is not None:
            plt.plot(alm_sp["r_norm"], alm_sp["fn"], "k--", label="ALM-LES")
        plt.xlabel("r/R")
        plt.ylabel("Fn (N/m)")
        plt.grid(True, alpha=0.3)
        plt.legend(title=study_param_meta["legend_title"], fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / "fn_by_cutoff.png", dpi=200)
        plt.close()

    if "ft" in args.plots:
        plt.figure(figsize=(6.5, 4.5))
        for r in spanwise:
            plt.plot(r["r_norm"], r["ft"], ".-", label=r["param_label"])
        if alm_sp is not None:
            plt.plot(alm_sp["r_norm"], alm_sp["ft"], "k--", label="ALM-LES")
        plt.xlabel("r/R")
        plt.ylabel("Ft (N/m)")
        plt.grid(True, alpha=0.3)
        plt.legend(title=study_param_meta["legend_title"], fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / "ft_by_cutoff.png", dpi=200)
        plt.close()

    print("Wrote plots to", out_dir)
    print(f"ALM baseline: Cp={cp_alm_mean:.4f}, Ct={ct_alm_mean:.4f}")


if __name__ == "__main__":
    main()
