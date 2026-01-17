import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Ensure we can import fvw
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from fvw.io import WakeReader

def main():
    parser = argparse.ArgumentParser(description="FVW Spanwise Analysis Tool")
    parser.add_argument("file", help="Path to wake.h5")
    parser.add_argument("--step", type=int, default=-1, help="Timestep to analyze (default: last step)")
    parser.add_argument("--blade", type=int, default=0, help="Blade index to analyze (default: 0)")
    args = parser.parse_args()

    h5_path = args.file
    if not os.path.exists(h5_path):
        print(f"Error: File {h5_path} not found.")
        return

    print(f"==================================================")
    print(f"   FVW Spanwise Analysis")
    print(f"   File: {h5_path}")
    print(f"==================================================")

    try:
        reader = WakeReader(h5_path)
        config = reader.read_config()
        timesteps = reader.list_timesteps()
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if not timesteps:
        print("No timesteps found in file.")
        return

    # Select timestep
    t_idx = args.step
    if t_idx == -1:
        # User wants last step.
        # Expert Tip: Avoid the very last step due to boundary artifacts. Use second to last.
        if len(timesteps) > 2:
            t_idx = timesteps[-2]
            print(f"Note: Auto-selecting second-to-last timestep {t_idx} to avoid boundary artifacts.")
        else:
            t_idx = timesteps[-1]
    else:
        if t_idx not in timesteps:
            print(f"Error: Timestep {t_idx} not found available steps: {timesteps[0]}..{timesteps[-1]}")
            return

    print(f"Analyzing Timestep: {t_idx}")

    # Read Data
    blade_data = reader.read_timestep_data(t_idx)
    if not blade_data or args.blade >= len(blade_data):
        print(f"Error: Invalid blade index or no data.")
        return

    bd = blade_data[args.blade] # Dict with cl, cd, aoa, rel_vel
    
    # Geometry
    r = config['r_center']
    R = config.get('r_tip', r[-1])
    r_R = r / R

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Subplot 1: Cl & Cd
    color = 'tab:blue'
    ax1.set_ylabel('Lift Coefficient ($C_l$)', color=color, fontsize=12)
    ax1.plot(r_R, bd['cl'], 'o-', color=color, label='$C_l$')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f"Spanwise Distribution at t={t_idx} (Blade {args.blade})")

    ax1b = ax1.twinx()
    color = 'tab:red'
    ax1b.set_ylabel('Drag Coefficient ($C_d$)', color=color, fontsize=12)
    ax1b.plot(r_R, bd['cd'], 's--', color=color, label='$C_d$')
    ax1b.tick_params(axis='y', labelcolor=color)

    # Subplot 2: AoA
    color = 'tab:green'
    ax2.set_xlabel('Normalized Radius ($r/R$)', fontsize=12)
    ax2.set_ylabel('Angle of Attack ($deg$)', color=color, fontsize=12)
    ax2.plot(r_R, bd['aoa'], '^-', color=color, label='$\\alpha$')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.grid(True, alpha=0.3)
    
    # Add Indication of typical Stall region (simplified)
    # ax2.axhline(y=15, color='gray', linestyle=':', label='Approx Stall')

    plt.tight_layout()
    
    # Save to the same directory as the input H5 file
    out_dir = os.path.dirname(h5_path)
    out_png = os.path.join(out_dir, f"spanwise_t{t_idx}.png")
    
    plt.savefig(out_png, dpi=150)
    print(f"\n[Success] Plot saved to: {os.path.abspath(out_png)}")

if __name__ == "__main__":
    main()
