import os
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import argparse

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description='Analyze ALM-LES Simulation Results (Ensight Format)')
parser.add_argument('case_dir', type=str, help='Path to the ALM-LES case directory')
parser.add_argument('--tasks', nargs='+', default=None, 
                    choices=['all', 'vor_z', 'vor_mag', 'inst_def', 'avg_def'],
                    help='Analysis tasks to run (Overridden by default behavior if not set)')
parser.add_argument('--avg', action='store_true', help='Run time-averaged analysis only (Default behavior runs instantaneous if this is not set)')

# Plot Limits
parser.add_argument('--w_min', type=float, default=-300, help='Vorticity Min')
parser.add_argument('--w_max', type=float, default=300, help='Vorticity Max')
parser.add_argument('--u_min', type=float, default=None, help='Velocity Min')
parser.add_argument('--u_max', type=float, default=None, help='Velocity Max')

# View Limits (in Diameter D)
parser.add_argument('--xlim', type=float, nargs=2, default=[-1.0, 5.0], help='X-axis limits in D (e.g. -1 5)')
parser.add_argument('--ylim', type=float, nargs=2, default=[-1.0, 1.0], help='Y-axis limits in D (e.g. -1 1)')

# Physics Parameters
parser.add_argument('--D', type=float, default=0.894, help='Rotor Diameter (m). Default: 0.894 (Scale Model)')
parser.add_argument('--U_inf', type=float, default=10.0, help='Free stream velocity (m/s). Default: 10.0')

args = parser.parse_args()

# --- Apply Default Task Logic ---
if args.tasks is None:
    if args.avg:
        args.tasks = ['avg_def']
    else:
        args.tasks = ['vor_z', 'vor_mag', 'inst_def']
    print(f"Running Tasks: {args.tasks}")

# Constants
D = args.D
U_INF = args.U_inf
CASE_DIR = args.case_dir
OUTPUT_DIR = os.path.join(CASE_DIR, 'analysis_output_les')
# Path to surfacesInst/horizontal/horizontal.case
ENSIGHT_PATH = os.path.join(CASE_DIR, 'postProcessing', 'surfacesInst', 'horizontal', 'horizontal.case')

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Analyzing Case: {CASE_DIR}")
print(f"Ensight Path: {ENSIGHT_PATH}")
print(f"Output Dir: {OUTPUT_DIR}")

if not os.path.exists(ENSIGHT_PATH):
    raise RuntimeError(f"Ensight file not found: {ENSIGHT_PATH}")

# Load Dataset
print("Loading Ensight dataset...")
reader = pv.get_reader(ENSIGHT_PATH)
# Reader time values
time_values = reader.time_values
print(f"Found {len(time_values)} timesteps. Range: {time_values[0]} to {time_values[-1]}")

# Helper to plot scalar
def plot_scalar(grid, scalar_name, filename, title, label, time_val=None, vmin=None, vmax=None, cmap='viridis'):
    points = grid.points
    # Rotate if necessary? 
    # Usually OpenFOAM: X=streamwise, Y=lateral, Z=vertical
    # For horizontal slice, Z is constant. Plot X vs Y.
    x = points[:, 0]
    y = points[:, 1]
    
    data = grid[scalar_name]
    
    plt.figure(figsize=(12, 4))
    
    if vmin is not None and vmax is not None:
        levels = np.linspace(vmin, vmax, 100)
    else:
        levels = 100
        
    plt.tricontourf(x/D, y/D, data, levels=levels, cmap=cmap, extend='both')
    
    cbar_ticks = None
    if vmin is not None and vmax is not None:
         cbar_ticks = np.linspace(vmin, vmax, 5)

    plt.colorbar(label=label, ticks=cbar_ticks)
    
    if time_val is not None:
         t_str = f"{time_val:.4f}"
    else:
         t_val = grid.field_data['TimeValue'][0] if 'TimeValue' in grid.field_data else 0.0
         t_str = f"{t_val:.4f}"
    
    plt.title(f'{title}\nt={t_str}s')
    plt.xlabel('x/D [-]')
    plt.ylabel('y/D [-]')
    plt.axis('scaled')
    
    if args.xlim is not None: plt.xlim(args.xlim)
    if args.ylim is not None: plt.ylim(args.ylim)
        
    plt.tight_layout()
    # Save
    out_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved {filename}")

def analyze_instantaneous():
    # Last timestep
    reader.set_active_time_value(time_values[-1])
    grid = reader.read()
    
    # It might be a MultiBlock, get the first block
    if isinstance(grid, pv.MultiBlock):
        grid = grid[0]
        
    print(f"Processing Instantaneous t={time_values[-1]}...")
    
    tasks = args.tasks
    if 'all' in tasks:
        tasks = ['vor_z', 'vor_mag', 'inst_def']
        
    # check compute derivatives
    if any(t in tasks for t in ['vor_z', 'vor_mag']):
        grid = grid.compute_derivative(scalars="U", vorticity=True)
        vorticity = grid['vorticity']
        grid['omega_z'] = vorticity[:, 2]
        grid['omega_mag'] = np.linalg.norm(vorticity, axis=1)
        
    if 'vor_z' in tasks:
        plot_scalar(grid, 'omega_z', 'les_inst_vor_z.png', 'Instantaneous Vertical Vorticity', 'Vertical Vorticity $\omega_z$ [1/s]', 
                   time_val=time_values[-1], vmin=args.w_min, vmax=args.w_max, cmap='RdBu_r')
                   
    if 'vor_mag' in tasks:
        plot_scalar(grid, 'omega_mag', 'les_inst_vor_mag.png', 'Instantaneous Vorticity Magnitude', 'Vorticity Magnitude $|\omega|$ [1/s]',
                   time_val=time_values[-1], vmin=0, vmax=args.w_max, cmap='inferno')
                   
    if 'inst_def' in tasks:
        U = grid['U']
        Ux = U[:, 0]
        deficit = (U_INF - Ux) / U_INF
        grid['deficit'] = deficit
        
        # Default 0-0.8 for deficit
        u_min = args.u_min if args.u_min is not None else 0.0
        u_max = args.u_max if args.u_max is not None else 1.0
        
        plot_scalar(grid, 'deficit', 'les_inst_deficit.png', 'Instantaneous Velocity Deficit', 'Instantaneous Velocity Deficit [-]',
                   time_val=time_values[-1], vmin=u_min, vmax=u_max, cmap='viridis')

def analyze_time_average():
    tasks = args.tasks
    if not ('all' in tasks or 'avg_def' in tasks):
        return
        
    print("Computing Time Average (Last 10% of timesteps)...")
    # Take last 10% or so
    n_steps = len(time_values)
    start_idx = max(0, n_steps - 50) # Arbitrary last 50 steps? Or last 10s?
    # Let's use last 100 if available, or all if less
    start_idx = max(0, n_steps - 100)
    
    subset_times = time_values[start_idx:]
    print(f"Averaging {len(subset_times)} steps from t={subset_times[0]} to {subset_times[-1]}")
    
    sum_Ux = None
    count = 0
    ref_grid = None
    
    for t in subset_times:
        reader.set_active_time_value(t)
        mesh = reader.read()
        if isinstance(mesh, pv.MultiBlock): mesh = mesh[0]
        
        Ux = mesh['U'][:, 0]
        
        if sum_Ux is None:
            sum_Ux = np.zeros_like(Ux)
            ref_grid = mesh
            
        sum_Ux += Ux
        count += 1
        
    avg_Ux = sum_Ux / count
    def_avg = (U_INF - avg_Ux) / U_INF
    ref_grid['def_avg'] = def_avg
    
    u_min = args.u_min if args.u_min is not None else 0.0
    u_max = args.u_max if args.u_max is not None else 1.0
    
    plot_scalar(ref_grid, 'def_avg', 'les_avg_deficit.png', 'Time-Averaged Velocity Deficit', 'Time-Averaged Velocity Deficit [-]',
               time_val=None, vmin=u_min, vmax=u_max, cmap='viridis')

if __name__ == "__main__":
    analyze_instantaneous()
    analyze_time_average()
    print("Done.")
