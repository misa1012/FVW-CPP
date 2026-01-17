import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# Ensure we can import fvw
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from fvw.io import WakeReader
from fvw.calc import PerformanceCalculator

def main():
    # Path to your H5 file - adjust as needed or pass as arg
    if len(sys.argv) > 1:
        h5_path = sys.argv[1]
    else:
        # Default fallback for testing
        h5_path = "../build/results/default_baseline/wake.h5"

    print(f"Analyzing: {h5_path}")
    
    try:
        reader = WakeReader(h5_path)
    except FileNotFoundError:
        print("File not found. Please run the simulation first.")
        return

    calc = PerformanceCalculator(reader)
    
    # Calculate stats. If simulation is short, start from 0.0
    res = calc.compute_time_series(start_rev=0.5)
    
    if res is None or len(res['time']) == 0:
        print("No valid data found or simulation too short.")
        return

    # Print mean stats
    mean_cp = np.mean(res['CP'])
    std_cp = np.std(res['CP'])
    mean_ct = np.mean(res['CT'])
    std_ct = np.std(res['CT'])
    
    print(f"\nResults ({len(res['time'])} steps examined):")
    print(f"Mean CP: {mean_cp:.4f} +/- {std_cp:.4f}")
    print(f"Mean CT: {mean_ct:.4f} +/- {std_ct:.4f}")
    
    # Simple Plot
    try:
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        color = 'tab:blue'
        ax1.set_xlabel('Rotor Revolutions')
        ax1.set_ylabel('Cp', color=color)
        ax1.plot(res['rev'], res['CP'], color=color, label='Cp')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, alpha=0.3)
        
        ax2 = ax1.twinx() 
        color = 'tab:orange'
        ax2.set_ylabel('Ct', color=color)  
        ax2.plot(res['rev'], res['CT'], color=color, label='Ct')
        ax2.tick_params(axis='y', labelcolor=color)
        
        plt.title(f"Aerodynamic Performance: {os.path.basename(os.path.dirname(h5_path))}")
        plt.tight_layout()
        
        out_png = "cp_ct_plot.png"
        plt.savefig(out_png)
        print(f"\nPlot saved to {out_png}")
        
    except ImportError:
        print("Matplotlib not installed, skipping plot.")

if __name__ == "__main__":
    main()
