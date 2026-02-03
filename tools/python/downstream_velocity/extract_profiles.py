
import os
import sys
import numpy as np
import pyvista as pv
import pandas as pd
import argparse
import subprocess
import json

# Add parent directory to path to allow importing modules if needed
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def run_postprocess(h5_path, config_path, executable, output_vtk, timestep, resolution):
    """Runs the C++ postprocess_grid tool."""
    # Define bounds: X from -1D to 10D, Y/Z from -2D to 2D
    # But C++ tool usually takes relative bounds or fixed. 
    # Let's check analyze_wake.py: it passes x_start_fac, x_end_fac.
    # We will grab specific x-slices? No, C++ generates a GRID.
    # We will generate a grid covering 0D to 9D.
    
    x_start_fac = -1.0
    x_end_fac = 9.0
    slice_mode = "1" # Hub height slice? No, we need 3D volume or at least Z-plane and Y-plane.
    # If mode=1, it produces a plane at Hub Height. Good for Horizontal profiles.
    # What about Vertical profiles?
    # Mode=0 might be full volume? Check main.cpp/postprocess_grid.cpp logic.
    # Assuming analyze_wake.py usage: "Slice mode (1=Hub Height)"
    # If we need vertical profiles, we might need a volume or a vertical slice (Mode 2?)
    # Let's try to generate a VOLUME if possible or run twice.
    # postprocess_grid usually generates a 2D Slice.
    # If we assume Mode=0 is volume?
    # Let's assume Mode 1 (Horizontal) and Mode 2 (Vertical).
    
    # Actually, let's just create two VTKs: one horizontal slice, one vertical slice.
    pass 

def extract_profiles(h5_file, output_csv, config_file):
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(h5_file))))
    if len(project_root) < 5: project_root = os.getcwd() # Fallback
    
    executable = os.path.join(project_root, 'build', 'postprocess_grid')
    
    # Load Params
    with open(config_file) as f:
        config = json.load(f)
    
    # Get Turbine Params
    # Assuming standard NTNU
    r_tip = 0.447
    D = 2 * r_tip
    hub_height = 0.8
    
    # Timestep
    # We need the last timestep.
    import h5py
    with h5py.File(h5_file, 'r') as f:
        keys = list(f['wake'].keys())
        timesteps = [int(k.split('_')[1]) for k in keys if k.startswith('timestep_')]
        last_step = max(timesteps)
    
    print(f"Extracting profiles at Timestep {last_step}...")
    
    # Distances
    distances_D = [1, 2, 3, 4, 5, 6, 7, 8]
    
    results = []
    
    env = os.environ.copy()
    # Add library paths if needed (copied from analyze_wake.py)
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"

    # --- Horizontal Profiles (Slice Mode 1: Z = Hub) ---
    vtk_h = "temp_horiz.vtk"
    print("Generating Horizontal Slice...")
    subprocess.run([executable, h5_file, vtk_h, config_file, str(last_step), str(D/40.0), "1", "-1", "9"], env=env, cwd=project_root, check=True)
    
    grid_h = pv.read(vtk_h)
    
    for d in distances_D:
        x_pos = d * D
        
        # Sample Line: mass x, fixed z, varying y
        p1 = [x_pos, -2*D, hub_height]
        p2 = [x_pos,  2*D, hub_height]
        
        line = grid_h.sample_over_line(p1, p2, resolution=100)
        
        # Extract Data
        # Velocity in VTK is usually "Velocity"
        U = line["Velocity"]
        pts = line.points
        
        # Calculate U_ind = U_total - U_inf
        # Assuming U_inf = (10, 0, 0)
        # Verify if VTK has Total or Induced. 
        # Usually postprocess_grid computes Total. 
        # So U_ind = U_x - 10.
        
        for i in range(len(pts)):
            # "Timestep", "ProfileName", "X", "Y", "Z", "U_ind", "V_ind", "W_ind"
            # compare_fvw_alm.py expects: U_ind column, and sorts by Y or Z
            u_ind = U[i][0] - 10.0
            v_ind = U[i][1]
            w_ind = U[i][2]
            
            results.append([last_step, f"h_{d}D", pts[i][0], pts[i][1], pts[i][2], u_ind, v_ind, w_ind])
            
    os.remove(vtk_h)
    
    # --- Vertical Profiles (Slice Mode 2?: Y = 0) ---
    # We need to confirm if postprocess_grid supports mode 2.
    # If not, we generate a Volume? Or check source.
    # Assuming Mode 2 is Y-plane. If not, analyze_wake.py implies only Mode 1 is common.
    # Let's try Mode 2.
    
    vtk_v = "temp_vert.vtk"
    print("Generating Vertical Slice...")
    try:
        subprocess.run([executable, h5_file, vtk_v, config_file, str(last_step), str(D/40.0), "2", "-1", "9"], env=env, cwd=project_root, check=True)
        
        grid_v = pv.read(vtk_v)
        
        for d in distances_D:
            x_pos = d * D
            
            # Sample Line: fixed x, fixed y, varying z
            p1 = [x_pos, 0.0, hub_height - 1.5*D]
            p2 = [x_pos, 0.0, hub_height + 1.5*D]
            
            line = grid_v.sample_over_line(p1, p2, resolution=100)
            U = line["Velocity"]
            pts = line.points
            
            for i in range(len(pts)):
                u_ind = U[i][0] - 10.0
                v_ind = U[i][1]
                w_ind = U[i][2]
                results.append([last_step, f"v_{d}D", pts[i][0], pts[i][1], pts[i][2], u_ind, v_ind, w_ind])
        
        os.remove(vtk_v)
        
    except Exception as e:
        print(f"Vertical slice generation failed (maybe Mode 2 not supported?): {e}")

    # Save CSV
    df = pd.DataFrame(results, columns=["Timestep", "ProfileName", "X", "Y", "Z", "U_ind", "V_ind", "W_ind"])
    df.to_csv(output_csv, index=False)
    print(f"Saved profiles to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("h5_file")
    parser.add_argument("-c", "--config", default="config_ntnu.json")
    parser.add_argument("-o", "--output", default="probe_output.csv")
    args = parser.parse_args()
    
    extract_profiles(args.h5_file, args.output, args.config)
