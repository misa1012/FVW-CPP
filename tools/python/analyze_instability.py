
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import json

def analyze_instability(h5_file, config_file="config_ntnu.json"):
    print(f"Analyzing {h5_file}...")
    
    # 1. Load Config to get nSegments
    try:
        with open(config_file, 'r') as f:
            cf = json.load(f)
            # Handle nested structure or flat
            if "turbine" in cf:
                n_segments = cf["turbine"].get("nSegments", 20)
            else:
                n_segments = cf.get("nSegments", 20)
    except:
        n_segments = 18 # Default fallback for NTNU
        print("Warning: Could not load config, defaulting nSegments=18")

    n_trail = n_segments + 1
    stride = n_trail
    print(f"Assuming Stride (nTrail) = {n_trail}")

    try:
        f = h5py.File(h5_file, 'r')
    except Exception as e:
        print(f"Error opening H5 file: {e}")
        return

    timesteps = sorted([int(k.split('_')[1]) for k in f['wake'].keys()])
    if not timesteps:
        return 
    last_t = timesteps[-1]
    print(f"Using Timestep: {last_t}")
    
    grp = f[f'wake/timestep_{last_t}']
    n_blades = len(grp.keys())
    
    tip_filaments = [] 
    hub_center = np.array([0, 0, 90.0]) # Or detect
    
    for b in range(n_blades):
        b_grp = grp[f'blade_{b}']
        nodes = b_grp['nodes'][:] # [N, 6]
        pos = nodes[:, 0:3]
        
        # Identify which "column" is the tip
        # We have 'stride' columns.
        # We check the average radius of each column to find the tip.
        
        # The column with max average radius might vary due to roll-up. 
        # For analysis, we SHOULD usually track the geometric tip node (the last one).
        # Geometric Tip = stride - 1 (since 0-indexed)
        tip_idx = stride - 1
        print(f"Blade {b}: Forcing Tip Index = {tip_idx} (Geometric Tip)")
        
        tip_nodes = pos[tip_idx::stride]
        
        # Sort by X? 
        # In FVW, nodes are Lagrangian. 
        # Typically the order in list is Time 0 -> Time Current.
        # So X should be decreasing (Oldest far downstream, Newest near rotor)?
        # Or Increasing?
        # Let's just sort by X to be safe for plotting lines.
        idx = np.argsort(tip_nodes[:,0])
        tip_nodes = tip_nodes[idx]
        
        tip_filaments.append(tip_nodes)

    f.close()
    
    # --- Plotting ---
    # Derive output path from input file
    # Input: results/CASE/wake.h5
    # Output: results/CASE/post_processing/instability_analysis
    case_dir = os.path.dirname(h5_file)
    output_dir = os.path.join(case_dir, "post_processing", "instability_analysis")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Outputting to: {output_dir}")
    R = 0.45 

    # 1. Radial Evolution
    plt.figure(figsize=(10, 6))
    for b, nodes in enumerate(tip_filaments):
        x = nodes[:, 0]
        r = np.sqrt((nodes[:,1] - hub_center[1])**2 + (nodes[:,2] - hub_center[2])**2)
        
        # Filter numerical ejections for clearer plot?
        # Only plot where r < 3*R
        mask = r < 3.0 * R
        plt.plot(x[mask]/R, r[mask]/R, '.', markersize=1, label=f'Blade {b}')
        
    plt.xlabel('x/R')
    plt.ylabel('r/R')
    plt.title('Tip Vortex Radial Position (No Filter)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/radial_instability.png")
    
    # 2. Pairing Distance
    if n_blades >= 2:
        b0 = tip_filaments[0]
        b1 = tip_filaments[1]
        
        # Filter ejections
        mask0 = (np.sqrt(b0[:,1]**2 + (b0[:,2]-90)**2) < 3.0*R)
        b0 = b0[mask0]
        
        mask1 = (np.sqrt(b1[:,1]**2 + (b1[:,2]-90)**2) < 3.0*R)
        b1 = b1[mask1]

        # Common X
        min_x = max(np.min(b0[:,0]), np.min(b1[:,0]))
        max_x = min(np.max(b0[:,0]), np.max(b1[:,0]))
        
        x_eval = np.linspace(min_x, max_x, 1000)
        
        y0 = np.interp(x_eval, b0[:,0], b0[:,1])
        z0 = np.interp(x_eval, b0[:,0], b0[:,2])
        y1 = np.interp(x_eval, b1[:,0], b1[:,1])
        z1 = np.interp(x_eval, b1[:,0], b1[:,2])
        
        dist = np.sqrt((y1-y0)**2 + (z1-z0)**2)
        
        plt.figure(figsize=(10, 6))
        plt.plot(x_eval/R, dist/R, 'k-', lw=1)
        plt.xlabel('x/R')
        plt.ylabel('Separation Distance / R')
        plt.title('Vortex Pairing Distance')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"{output_dir}/pairing_distance.png")

if __name__ == "__main__":
    analyze_instability("results/NTNU_Baseline/wake.h5")
