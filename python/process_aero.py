import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# Ensure we can import fvw
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from fvw.io import WakeReader
from fvw.calc import PerformanceCalculator

def main():
    """
    FVW Aerodynamic Analysis Tool
    
    Usage: python process_aero.py [path_to_h5_file]
    """
    
    # 1. Parse Arguments
    if len(sys.argv) > 1:
        h5_path = sys.argv[1]
    else:
        # Smart default: look for the most recent result
        default_dir = "../results/case_1"
        h5_path = os.path.join(default_dir, "wake.h5")
        if not os.path.exists(h5_path):
             print(f"Usage: python process_aero.py <path/to/wake.h5>")
             print(f"Error: Default file {h5_path} not found.")
             return

    print(f"==================================================")
    print(f"   FVW Aerodynamic Analysis")
    print(f"   File: {h5_path}")
    print(f"==================================================")
    
    # 2. Load Data
    try:
        reader = WakeReader(h5_path)
    except FileNotFoundError:
        print("Error: File not found.")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    calc = PerformanceCalculator(reader)
    config = reader.read_config()
    
    print(f"Sim Configuration:")
    print(f"  - Wind Speed: {config.get('Uinf', 'N/A')} m/s")
    print(f"  - RPM: {config.get('omega', 0) * 9.5493:.2f} rpm")
    print(f"--------------------------------------------------")

    # 3. Compute (Start from 0.5 revs to skip initial shock)
    try:
        res = calc.compute_time_series(start_rev=0.5)
    except Exception as e:
        print(f"Calculation failed: {e}")
        return
    
    if res is None or len(res['time']) == 0:
        print("No valid data found or simulation too short.")
        return

    # 4. Filter Artifacts (Remove last 2 steps)
    # The last step often has boundary artifacts due to backward difference.
    # The user requested removing the last 2 steps to be safe.
    cutoff = 2
    if len(res['time']) > cutoff:
        print(f"Filtering: Removing last {cutoff} timesteps (boundary artifacts).")
        for key in res:
            if isinstance(res[key], (np.ndarray, list)):
                res[key] = res[key][:-cutoff]
    else:
        print("Warning: Data series too short to filter artifacts.")

    # 5. Statistics
    mean_cp = np.mean(res['CP'])
    std_cp = np.std(res['CP'])
    mean_ct = np.mean(res['CT'])
    std_ct = np.std(res['CT'])
    mean_power = np.mean(res['Power']) / 1e6 # MW
    mean_thrust = np.mean(res['Thrust']) / 1e3 # kN
    
    print(f"\nResults ({len(res['time'])} steps):")
    print(f"  Mean Cp:     {mean_cp:.4f} +/- {std_cp:.4f}")
    print(f"  Mean Ct:     {mean_ct:.4f} +/- {std_ct:.4f}")
    print(f"  Mean Power:  {mean_power:.3f} MW")
    print(f"  Mean Thrust: {mean_thrust:.3f} kN")
    
    # 6. Plotting
    try:
        # Figure 1: Cp & Ct
        fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
        
        # Plot Cp
        color = 'tab:blue'
        ax1.set_ylabel('Power Coefficient ($C_P$)', color=color, fontsize=12)
        ax1.plot(res['rev'], res['CP'], color=color, linewidth=2, label='$C_P$')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, alpha=0.3)
        ax1.set_title(f"Aerodynamic Coefficients: {os.path.basename(os.path.dirname(h5_path))}")
        
        # Twin axis for Ct
        ax2 = ax1.twinx() 
        color = 'tab:orange'
        ax2.set_ylabel('Thrust Coefficient ($C_T$)', color=color, fontsize=12)  
        ax2.plot(res['rev'], res['CT'], color=color, linewidth=2, linestyle='--', label='$C_T$')
        ax2.tick_params(axis='y', labelcolor=color)
        
        # Plot Power & Thrust (Absolute)
        color = 'tab:green'
        ax3.set_xlabel('Rotor Revolutions', fontsize=12)
        ax3.set_ylabel('Power [MW]', color=color, fontsize=12)
        ax3.plot(res['rev'], res['Power']/1e6, color=color, linewidth=2, label='Power')
        ax3.tick_params(axis='y', labelcolor=color)
        ax3.grid(True, alpha=0.3)
        
        ax4 = ax3.twinx()
        color = 'tab:red'
        ax4.set_ylabel('Thrust [kN]', color=color, fontsize=12)
        ax4.plot(res['rev'], res['Thrust']/1e3, color=color, linewidth=2, linestyle='--', label='Thrust')
        ax4.tick_params(axis='y', labelcolor=color)

        plt.tight_layout()
        
        # Save to the same directory as the input H5 file
        out_dir = os.path.dirname(h5_path)
        out_png = os.path.join(out_dir, "aero_performance.png")
        plt.savefig(out_png, dpi=150)
        print(f"\n[Success] Plot saved to: {os.path.abspath(out_png)}")
        
    except ImportError:
        print("\n[Warning] Matplotlib not installed, skipping plot.")

if __name__ == "__main__":
    main()
