import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

# Add parent directory to sys.path to access utils if needed
# script_dir = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(os.path.join(script_dir, '..', '..'))

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize ALM-LES Wake Fields from VTK slices.")
    parser.add_argument("--alm_dir", required=True, help="Path to ALM-LES case directory")
    parser.add_argument("--time", default="1.872", help="Time step to analyze (default: 1.872)")
    parser.add_argument("--u_inf", type=float, default=10.0, help="Freesteam velocity (default: 10.0)")
    parser.add_argument("--diameter", type=float, default=0.894, help="Rotor Diameter D (default: 0.894)")
    parser.add_argument("--hub_height", type=float, default=0.8, help="Hub Height (default: 0.8)")
    parser.add_argument("--output_dir", default="results/alm_wake_plots", help="Output directory for plots")
    return parser.parse_args()

def plot_yz_deficit(grid, label, args):
    """
    Plots YZ velocity deficit contour from a PyVista grid (slice).
    """
    # Grid points
    points = grid.points # (N, 3)
    # Velocity field 'U' -> (N, 3)
    U = grid.point_data['U']
    Ux = U[:, 0] # Streamwise component
    
    # Calculate Deficit
    deficit = 1.0 - (Ux / args.u_inf)
    
    # Coordinates
    y = points[:, 1]
    z = points[:, 2]
    
    # Normalize coordinates
    y_norm = y / args.diameter
    z_norm = (z - args.hub_height) / args.diameter
    
    print(f"    Data Range -> y_norm: [{y_norm.min():.3f}, {y_norm.max():.3f}], z_norm: [{z_norm.min():.3f}, {z_norm.max():.3f}]")
    
    # Create triangulation for plotting
    # Since VTK data is unstructured/polydata, simple reshape might not work directly 
    # unless it's a structured grid. The 'surfaces' function usually outputs PolyData (triangles/polys).
    # Tricontourf is robust for this.
    
    plt.figure(figsize=(8, 6))
    
    # contourf
    # Match calculate_rews.py style: viridis, 0 to 1.0
    levels = np.linspace(0, 1.0, 50) 
    tc = plt.tricontourf(y_norm, z_norm, deficit, levels=levels, cmap='viridis', extend='both')
    cbar = plt.colorbar(tc, label=r"Velocity Deficit ($1 - U_x/U_\infty$)")
    
    # Add Rotor Circle
    # R normalized by D is 0.5. Center is (0, 0) in normalized coords.
    theta = np.linspace(0, 2*np.pi, 200)
    r_circle = 0.5
    xc = r_circle * np.cos(theta)
    zc = r_circle * np.sin(theta)
    plt.plot(xc, zc, 'r--', linewidth=2, label='Rotor Tip')
    
    # Labels and Limits
    plt.xlabel(r"$y/D$")
    plt.ylabel(r"$(z - z_{hub})/D$")
    plt.title(f"ALM-LES Wake Deficit at {label} (t={args.time}s)")
    plt.xlim([-1.0, 1.0])
    plt.ylim([-0.78, 1.0]) # Crop bottom empty region (mesh starts at ~-0.79D)
    plt.gca().set_aspect('equal', adjustable='box') # Enforce square aspect ratio without changing limits
    plt.grid(True, alpha=0.3)
    
    # Save
    out_path = os.path.join(args.output_dir, f"alm_deficit_{label}.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {out_path}")
    plt.close()

def main():
    args = parse_args()
    
    # VTK Directory Detection
    # Convention: postProcessing/surfaces_yz/<time>/<name>.vtk OR postProcessing/surfaces/...
    # Search for logical candidates
    candidates = ["surfaces_yz", "surfaces", "surfaces_yz_flat", "surfaces_yz_corrected"]
    vtk_dir = None
    
    for c in candidates:
        d = os.path.join(args.alm_dir, "postProcessing", c, args.time)
        if os.path.exists(d):
            vtk_dir = d
            print(f"Found VTK directory: {vtk_dir}")
            break
            
    if not vtk_dir:
        print(f"Error: Could not find VTK directory for time {args.time} in {args.alm_dir}/postProcessing/")
        print(f"Checked: {candidates}")
        print("Did you run 'postProcess -func ...'?")
        sys.exit(1)
        
    os.makedirs(args.output_dir, exist_ok=True)
    
    # List all VTK/VTP files
    # Expected names: yz_1D.vtp, etc.
    import glob
    files = glob.glob(os.path.join(vtk_dir, "*_?D.vtk")) + glob.glob(os.path.join(vtk_dir, "*_?D.vtp"))
    
    if not files:
        print(f"No VTK/VTP files found in {vtk_dir}")
        sys.exit(1)
        
    print(f"Found {len(files)} slice files. Processing...")
    
    for fpath in files:
        name = os.path.basename(fpath).split('.')[0] # remove extension
        if name.endswith("_U"): name = name[:-2] # handle legacy _U naming if present
        
        print(f" - Processing {name} ({fpath})...")
        
        try:
            grid = pv.read(fpath)
            plot_yz_deficit(grid, name, args)
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    main()
