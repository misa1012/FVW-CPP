
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import argparse

def parse_xy_file(filepath):
    """Parses OpenFOAM .xy file."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 4: # y, x, z. u_x, u_y, u_z 
                # ALM sample output format typically: (x y z) (ux uy uz)
                # But sometimes just coord value and vector.
                # Let's check format. It seems to be .xy from 'sample' utility.
                # Usually coord, vector.
                pass
            data.append([float(x) for x in parts[0:4]]) # x, y, z, val (scalar) or vec
            # Wait, UMean.xy likely has (x y z) (ux uy uz)
            
    # Re-reading using numpy is safer if format is consistent
    # format: y U_x U_y U_z (if sample set is line)
    return np.loadtxt(filepath) 

def load_alm_data(alm_dir, time_dir="1.872"):
    """Loads ALM data from specific time directory."""
    path = os.path.join(alm_dir, "postProcessing/sample", time_dir)
    data = {}
    
    # Horizontal profiles
    # format: horizontal_1D_UMean.xy
    # Content typically: distance_along_line (x y z) (u v w) ??
    # Actually OpenFOAM sample sets usually output: x y z u v w for vector fields
    
    # Let's inspect one file first to be sure about column indices
    # We will assume columns: 0:x, 1:y, 2:z, 3:ux, 4:uy, 5:uz
    
    files = glob.glob(os.path.join(path, "horizontal_*_UMean.xy"))
    for f in files:
        fname = os.path.basename(f)
        try:
            # Extract x/D location from filename, e.g., horizontal_1D_UMean.xy -> 1D
            label = fname.split('_')[1] # 1D, 2D...
            
            arr = np.loadtxt(f)
            # Filter: we need normalized U/U_inf. U_inf = 10.0 (usually)
            # Store (y, u_x) for horizontal
            data[f"h_{label}"] = {
               "y": arr[:, 1], # Assuming y varies
               "u": arr[:, 3]  # Assuming ux is col 3
            }
        except Exception as e:
            print(f"Skipping {fname}: {e}")

    # Vertical profiles
    files = glob.glob(os.path.join(path, "vertical_*_UMean.xy"))
    for f in files:
        fname = os.path.basename(f)
        try:
            label = fname.split('_')[1]
            arr = np.loadtxt(f)
            data[f"v_{label}"] = {
               "z": arr[:, 2], # Assuming z varies
               "u": arr[:, 3]
            }
        except Exception as e:
            print(f"Skipping {fname}: {e}")
            
    return data

def load_fvw_data(csv_path, last_timestep=None):
    df = pd.read_csv(csv_path)
    if last_timestep:
        df = df[df['Timestep'] == last_timestep]
    else:
        last_timestep = df['Timestep'].max()
        df = df[df['Timestep'] == last_timestep]
    
    print(f"Loaded FVW data for timestep {last_timestep}")
    
    data = {}
    
    # Group by ProfileName
    # ProfileName format in FVW: h_1D, h_1_5D...
    # Need to match ALM format: h_1D matches horizontal_1D
    
    for name, group in df.groupby("ProfileName"):
        # Calculate U_total = U_free + U_ind
        # U_free = 10.0, 0, 0
        U_inf = 10.0
        
        # Sort by spatial coord
        if name.startswith("h_"):
            group = group.sort_values("Y")
            data[name] = {
                "y": group["Y"].values,
                "u": U_inf + group["U_ind"].values
            }
        elif name.startswith("v_"):
            group = group.sort_values("Z")
            data[name] = {
                "z": group["Z"].values,
                "u": U_inf + group["U_ind"].values
            }
            
    return data

def compare_and_plot(alm_data, fvw_data, output_dir=None):
    if output_dir is None:
        output_dir = "comparison_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    # Constants for Normalization
    HUB_ALM = 0.82 # Approx (from 0.894/2 + clearance?) or 0.8. NTNU usually 0.82.
    HUB_FVW = 90.0 # Detected from wake bounds
    R = 0.45       # Approx
    
    distances = ["1D", "2D", "3D", "4D", "5D", "6D", "7D", "8D"]
    
    for d in distances:
        # Horizontal
        alm_key = f"h_{d}"
        fvw_key = f"h_{d}"
        
        if alm_key in alm_data and fvw_key in fvw_data:
            plt.figure(figsize=(6, 4))
            
            # ALM
            plt.plot(alm_data[alm_key]["y"]/R, alm_data[alm_key]["u"]/10.0, 'k--', label="ALM-LES (Mean)")
            
            # FVW
            plt.plot(fvw_data[fvw_key]["y"]/R, fvw_data[fvw_key]["u"]/10.0, 'r-', label="FVW (Inst.)")
            
            plt.xlabel("y/R")
            plt.ylabel("U/U_inf")
            plt.title(f"Horizontal Velocity Profile at x={d}")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{output_dir}/compare_horizontal_{d}.png")
            plt.close()
            print(f"Saved {output_dir}/compare_horizontal_{d}.png")

        # Vertical
        alm_key = f"v_{d}"
        fvw_key = f"v_{d}"
        
        if alm_key in alm_data and fvw_key in fvw_data:
            plt.figure(figsize=(6, 4))
            
            # ALM - Normalize Z
            z_alm_norm = (alm_data[alm_key]["z"] - HUB_ALM) / R
            
            plt.plot(alm_data[alm_key]["u"]/10.0, z_alm_norm, 'k--', label="ALM-LES (Mean)")
            
            # FVW - Normalize Z
            # Note: FVW 'z' column from csv is absolute (around 90)
            z_fvw_norm = (fvw_data[fvw_key]["z"] - HUB_FVW) / R
            
            plt.plot(fvw_data[fvw_key]["u"]/10.0, z_fvw_norm, 'r-', label="FVW (Inst.)")
            
            plt.xlabel("U/U_inf")
            plt.ylabel("(z - Hub)/R")
            plt.title(f"Vertical Velocity Profile at x={d}")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{output_dir}/compare_vertical_{d}.png")
            plt.close()
            print(f"Saved {output_dir}/compare_vertical_{d}.png")

if __name__ == "__main__":
    import sys
    # Hardcoded paths data for now
    alm_path = "/data/marine-cfd/shug8104/ALM-LES/03_v2_Refined_NoNacelle_nx150_nt240"
    fvw_csv = "results/NTNU_Baseline/probe_output.csv"
    
    # Deduce output directory
    case_dir = os.path.dirname(fvw_csv)
    out_dir = os.path.join(case_dir, "post_processing", "comparison_plots")
    
    print("Loading ALM data...")
    alm = load_alm_data(alm_path)
    
    # Check if CSV exists
    if not os.path.exists(fvw_csv):
        print(f"Error: {fvw_csv} not found. Run post-processing first.")
        # Try finding partial? No, user deleted it.
    else:
        print("Loading FVW data...")
        fvw = load_fvw_data(fvw_csv)
        
        print(f"Generating plots to {out_dir}...")
        compare_and_plot(alm, fvw, output_dir=out_dir)
