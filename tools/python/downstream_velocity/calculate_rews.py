
import os
import sys
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import subprocess
import json

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

def run_postprocess(h5_path, config_path, executable, output_vtk, timestep, resolution, x_loc_D):
    """Runs postprocess_grid to generate a YZ slice at specific X location."""
    # We use Slice Mode 0 (Volume) but restrict X to a single plane
    
    # x_start = x_end = x_loc_D
    
    cmd = [
        executable,
        h5_path,
        output_vtk,
        config_path,
        str(timestep),
        str(resolution),
        "0", # Mode 0 = Volume (Standard 3D grid)
        str(x_loc_D), # x_start
        str(x_loc_D)  # x_end
    ]
    
    # Environment setup
    env = os.environ.copy()
    # Add necessary library paths (generic hpc setup)
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"

    # Assume PROJECT_ROOT is 3 levels up from this script (tools/python/downstream_velocity)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "../../.."))
    
    subprocess.run(cmd, env=env, cwd=project_root, check=True, stdout=subprocess.DEVNULL)

def analyze_rews(h5_file, config_file, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(h5_file), "post_processing", "rews_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup Paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "../../.."))
    executable = os.path.join(project_root, 'build', 'postprocess_grid')
    
    # Load Config
    with open(config_file) as f:
        config = json.load(f)
        
    r_tip = 0.447
    D = 2.0 * r_tip
    hub_height = 0.8
    U_inf = 10.0
    
    # Get Max Timestep
    import h5py
    with h5py.File(h5_file, 'r') as f:
        keys = list(f['wake'].keys())
        timesteps = [int(k.split('_')[1]) for k in keys if k.startswith('timestep_')]
        last_step = max(timesteps)
        
    print(f"Analyzing REWS at t={last_step} (dt step)...")
    
    distances = [1, 2, 3, 4, 5]
    
    metrics = []
    
    print(f"{'Dist(D)':<10} | {'REWS (m/s)':<15} | {'Deficit (-)':<15}")
    print("-" * 45)
    
    for d in distances:
        vtk_name = f"temp_slice_{d}D.vtk"
        vtk_path = os.path.join(output_dir, vtk_name)
        
        # 1. Generate VTK Slice
        # Resolution: use reasonably fine grid for integration
        res = D / 40.0 # ~0.02m
        run_postprocess(h5_file, config_file, executable, vtk_path, last_step, res, d)
        
        # 2. Load Data
        grid = pv.read(vtk_path)
        
        # Points: (N, 3)
        pts = grid.points
        # Velocity: (N, 3)
        vel = grid["Velocity"]
        ux = vel[:, 0]
        
        y = pts[:, 1]
        z = pts[:, 2]
        
        # 3. Filter for Rotor Area
        # Distance from hub center (y=0, z=hub_height)
        r_dist = np.sqrt(y**2 + (z - hub_height)**2)
        mask = r_dist <= r_tip
        
        # 4. Calculate REWS
        if np.any(mask):
            u_rotor = ux[mask]
            # REWS = cuberoot( mean( u^3 ) )
            u_cubed_mean = np.mean(u_rotor**3)
            rews = np.cbrt(u_cubed_mean)
            deficit = 1.0 - (rews / U_inf)
        else:
            rews = np.nan
            deficit = np.nan
            
        print(f"{d:<10} | {rews:<15.3f} | {deficit:<15.3f}")
        metrics.append({"Distance_D": d, "REWS": rews, "Deficit": deficit})
        
        # 5. Visualize YZ Deficit Map
        plt.figure(figsize=(6, 5))
        
        # Calculate local deficit
        local_deficit = 1.0 - (ux / U_inf)
        
        # Grid/Contour
        # Use tricontourf for unstructured points (pyvista reads as PointSet effectively here)
        levels = np.linspace(0, 1.0, 50)
        # Normalize by D
        tc = plt.tricontourf(y/D, (z - hub_height)/D, local_deficit, levels=levels, cmap="viridis", extend="both")
        plt.colorbar(tc, label=r"Velocity Deficit $1 - U_x/U_{\infty}$")
        
        # Draw Rotor Circle (Red Dashed)
        # Normalized coordinates: y/D, (z-H)/D
        # Circle at (0,0) with radius 0.5 (since R = 0.5*D)
        circle = plt.Circle((0, 0), 0.5, color='red', fill=False, linestyle='--', linewidth=2, label="Rotor Plane")
        plt.gca().add_patch(circle)
        
        plt.xlabel("y/D")
        plt.ylabel("(z - Hub)/D")
        plt.title(f"Wake Cross-section at {d}D (REWS Deficit={deficit:.3f})")
        plt.axis('scaled')
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)
        # plt.legend() # Legend for patch might act weird with contour, usually not needed if explicit
        
        plot_path = os.path.join(output_dir, f"wake_yz_{d}D.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        
        # Cleanup
        os.remove(vtk_path)
        
    # Save Metrics CSV
    df = pd.DataFrame(metrics)
    csv_path = os.path.join(output_dir, "rews_metrics.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved metrics to {csv_path}")
    print(f"Saved plots to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_file")
    parser.add_argument("-c", "--config", default="config_ntnu.json")
    parser.add_argument("-o", "--output_dir", default=None)
    args = parser.parse_args()
    
    analyze_rews(args.h5_file, args.config, args.output_dir)
