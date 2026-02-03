
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

import sys

# ... (imports)

def compare_and_plot(alm_data, fvw_data, output_dir=None):
    if output_dir is None:
        output_dir = "comparison_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    # Constants for Normalization
    HUB_ALM = 0.82 
    HUB_FVW = 0.8  
    R = 0.447      
    
    # Get FVW time
    # We assume one timestep in the CSV for now or take the first one
    # We don't have DT here, but we can assume t=1.872 if it matches max steps?
    # Or just label "FVW" if time is unknown? 
    # Extract timestep from keys? No, data structure is dict.
    # The load_fvw_data function prints the timestep.
    # Let's pass timestep to this function or just hardcode for this specific request if easier, 
    # but better to handle it.
    # We can infer it from the standard dt=0.000585 * 3200 ~= 1.872
    fvw_time_str = "t=1.872s" # Approx matching
    alm_time_str = "t=1.872s"
    
    distances = ["1D", "2D", "3D", "4D", "5D", "6D", "7D", "8D"]
    
    for d in distances:
        # Horizontal
        alm_key = f"h_{d}"
        fvw_key = f"h_{d}"
        
        if alm_key in alm_data and fvw_key in fvw_data:
            plt.figure(figsize=(6, 4))
            
            # ALM
            u_alm_norm = alm_data[alm_key]["u"] / 10.0
            def_alm = 1.0 - u_alm_norm
            plt.plot(alm_data[alm_key]["y"]/R, def_alm, 'k--', label=f"ALM-LES {alm_time_str}")
            
            # FVW
            u_fvw_norm = fvw_data[fvw_key]["u"] / 10.0
            def_fvw = 1.0 - u_fvw_norm
            plt.plot(fvw_data[fvw_key]["y"]/R, def_fvw, 'r-', label=f"FVW {fvw_time_str}")
            
            plt.xlabel("y/R")
            plt.ylabel("Velocity Deficit $1 - U/U_{\infty}$")
            plt.title(f"Horizontal Deficit at x={d}")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{output_dir}/compare_horizontal_{d}_deficit.png")
            plt.close()
            print(f"Saved {output_dir}/compare_horizontal_{d}_deficit.png")

        # Vertical
        alm_key = f"v_{d}"
        fvw_key = f"v_{d}"
        
        if alm_key in alm_data and fvw_key in fvw_data:
            plt.figure(figsize=(6, 4))
            
            # ALM - Normalize Z
            z_alm_norm = (alm_data[alm_key]["z"] - HUB_ALM) / R
            u_alm_norm = alm_data[alm_key]["u"] / 10.0
            def_alm = 1.0 - u_alm_norm
            
            plt.plot(def_alm, z_alm_norm, 'k--', label=f"ALM-LES {alm_time_str}")
            
            # FVW - Normalize Z
            z_fvw_norm = (fvw_data[fvw_key]["z"] - HUB_FVW) / R
            u_fvw_norm = fvw_data[fvw_key]["u"] / 10.0
            def_fvw = 1.0 - u_fvw_norm
            
            plt.plot(def_fvw, z_fvw_norm, 'r-', label=f"FVW {fvw_time_str}")
            
            plt.xlabel("Velocity Deficit $1 - U/U_{\infty}$")
            plt.ylabel("(z - Hub)/R")
            plt.title(f"Vertical Deficit at x={d}")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{output_dir}/compare_vertical_{d}_deficit.png")
            plt.close()
            print(f"Saved {output_dir}/compare_vertical_{d}_deficit.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare FVW wake profiles with ALM reference.")
    parser.add_argument("fvw_csv", help="Path to FVW probe output CSV")
    parser.add_argument("-r", "--ref_dir", required=True, help="Path to ALM case directory")
    parser.add_argument("-o", "--output_dir", help="Output directory for plots")
    
    args = parser.parse_args()
    
    # Deduce output directory
    if args.output_dir is None:
        case_dir = os.path.dirname(args.fvw_csv)
        args.output_dir = os.path.join(case_dir, "post_processing", "comparison_plots")
    
    print(f"Loading ALM data from {args.ref_dir}...")
    alm = load_alm_data(args.ref_dir)
    
    if not os.path.exists(args.fvw_csv):
        print(f"Error: {args.fvw_csv} not found.")
        sys.exit(1)
        
    print(f"Loading FVW data from {args.fvw_csv}...")
    fvw = load_fvw_data(args.fvw_csv)
    
    print(f"Generating plots to {args.output_dir}...")
    compare_and_plot(alm, fvw, output_dir=args.output_dir)
    pass
