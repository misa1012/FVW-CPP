
import numpy as np
import pandas as pd
import glob
import os

def load_xy(filepath):
    try:
        arr = np.loadtxt(filepath)
        # Assuming col 1 is Y (or Z), col 3 is U_x
        # ALM format usually: y/z x val or x y z u v w
        # Based on previous check: 
        # horizontal: col 1 is y, col 3 is u
        # vertical: col 2 is z, col 3 is u
        return arr
    except:
        return None

def analyze():
    # Constants
    U_inf = 10.0
    R = 0.45
    Hub_ALM = 0.82
    Hub_FVW = 90.0 # From previous checks
    
    alm_path = "/data/marine-cfd/shug8104/ALM-LES/03_v2_Refined_NoNacelle_nx150_nt240/postProcessing/sample/1.872"
    fvw_csv = "results/NTNU_Baseline/probe_output_partial.csv"
    
    print("--- Fluid Dynamics Analysis Metrics ---")
    
    # 1. ALM Data
    print("\n[ALM-LES Results (t=1.872s)]")
    for d in [1, 3, 5, 7]:
        # Horizontal
        f = f"{alm_path}/horizontal_{d}D_UMean.xy"
        if os.path.exists(f):
            data = load_xy(f)
            if data is not None:
                u = data[:, 3]
                u_min = np.min(u)
                deficit = 1.0 - (u_min / U_inf)
                print(f"  {d}D Horizontal: Max Deficit = {deficit:.3f} (U_min={u_min:.2f})")
    
    # 2. FVW Data
    print("\n[FVW Results (t=1.12s)]")
    try:
        df = pd.read_csv(fvw_csv)
        final_t = df['Timestep'].max()
        df = df[df['Timestep'] == final_t]
        
        for d in [1, 3, 5, 7]:
            # Look for profile name h_{d}D
            # Note: script used h_1D, h_2D... need to match generation logic
            # Postprocess logic was: h_1D, h_1_5D etc.
            name = f"h_{d}D"
            subset = df[df['ProfileName'] == name]
            if not subset.empty:
                # U_ind is in U_ind column. U_total = 10 + U_ind
                u_ind = subset['U_ind'].values
                u_total = U_inf + u_ind
                u_min = np.min(u_total)
                deficit = 1.0 - (u_min / U_inf)
                
                # Check Centerline (y~0)
                # Filter y close to 0
                center = subset[subset['Y'].abs() < 0.1]
                u_center = U_inf + center['U_ind'].mean() if not center.empty else u_min
                
                print(f"  {d}D Horizontal: Max Deficit = {deficit:.3f} (U_min={u_min:.2f})")
    except Exception as e:
        print(f"Error loading FVW: {e}")

if __name__ == "__main__":
    analyze()
