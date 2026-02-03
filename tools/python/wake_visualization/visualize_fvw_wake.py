
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import json
import argparse
import h5py

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description='Analyze Wake Simulation Results')
parser.add_argument('--tasks', nargs='+', default=None, 
                    choices=['all', 'vor_z', 'vor_mag', 'inst_def', 'avg_def'],
                    help='Analysis tasks to run (Overridden by default behavior if not set)')
parser.add_argument('--avg', action='store_true', help='Run time-averaged analysis only (Default behavior runs instantaneous if this is not set)')

# Path settings
parser.add_argument('--results_dir', type=str, default=None, help='Results directory (default: <root>/results/NTNU_Baseline)')
parser.add_argument('--config', type=str, default=None, help='Config file path (default: <root>/config_ntnu.json)')

# Time settings (Defaults are calculated from config if not provided)
parser.add_argument('--t_inst', type=int, default=None, help='Timestep for instantaneous vorticity plot (Default: Last Step)')
parser.add_argument('--t_start', type=int, default=None, help='Start timestep for averaging (Default: Last 10 intervals)')
parser.add_argument('--t_end', type=int, default=None, help='End timestep for averaging (Default: Last Step)')
parser.add_argument('--avg_step', type=int, default=10, help='Step size for averaging')

# Plot Limits (Default: None = Auto-scale based on data range)
parser.add_argument('--w_min', type=float, default=-300, help='Vorticity Min')
parser.add_argument('--w_max', type=float, default=300, help='Vorticity Max')
parser.add_argument('--u_min', type=float, default=None, help='Velocity Min')
parser.add_argument('--u_max', type=float, default=None, help='Velocity Max')

# View Limits (in Diameter D)
parser.add_argument('--xlim', type=float, nargs=2, default=[-1.0, 5.0], help='X-axis limits in D (e.g. -1 5)')
parser.add_argument('--ylim', type=float, nargs=2, default=[-1.0, 1.0], help='Y-axis limits in D (e.g. -1 1)')

# Resolution & Slice settings
parser.add_argument('--res_vort', type=float, default=None, help='Grid resolution for vorticity (m, overrides ppd)')
parser.add_argument('--res_avg', type=float, default=None, help='Grid resolution for averaging (m, overrides ppd)')
parser.add_argument('--ppd_vort', type=float, default=160.0, help='Points per Diameter for vorticity (Default: 80)')
parser.add_argument('--ppd_avg', type=float, default=80.0, help='Points per Diameter for averaging (Default: 80)')
parser.add_argument('--slice_mode', type=str, default="1", help='Slice mode (1=Hub Height)')

args = parser.parse_args()

# --- Apply Default Task Logic ---
if args.tasks is None:
    if args.avg:
        args.tasks = ['avg_def']
    else:
        args.tasks = ['vor_z', 'vor_mag', 'inst_def']
    print(f"Running Tasks: {args.tasks}")

# --- Dynamic Path Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# SCRIPT_DIR is tools/python/wake_visualization/
# Project Root is ../../../
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

# Defaults relative to project root
DEFAULT_RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'NTNU_0_25')
DEFAULT_CONFIG = os.path.join(PROJECT_ROOT, 'config_ntnu.json')

CASE_DIR = args.results_dir if args.results_dir else DEFAULT_RESULTS_DIR
CONFIG_PATH = args.config if args.config else DEFAULT_CONFIG

BUILD_DIR = os.path.join(PROJECT_ROOT, 'build')
EXECUTABLE = os.path.join(BUILD_DIR, 'postprocess_grid')
H5_PATH = os.path.join(CASE_DIR, 'wake.h5')
OUTPUT_DIR = os.path.join(CASE_DIR, 'analysis_output')

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Project Root: {PROJECT_ROOT}")
print(f"Results Dir: {CASE_DIR}")
print(f"Config: {CONFIG_PATH}")

# --- Auto-Detect Defaults from Config & HDF5 ---
try:
    with open(CONFIG_PATH, 'r') as f:
        config_data = json.load(f)
    
    # Sim Params (Draft max_step from config)
    sim_config = config_data.get('simulation', {})
    n_revs = sim_config.get('numRevolutions', 10.0)
    steps_per_rev = sim_config.get('stepsPerRevolution', 100)
    max_step = int(n_revs * steps_per_rev)
    
    # Check actual HDF5 file for max available timestep
    # This overrides the config-based guess which might be wrong (e.g. if sim crashed)
    if os.path.exists(H5_PATH):
        try:
            with h5py.File(H5_PATH, 'r') as f:
                # Check root or 'wake' group
                target_group = f
                if 'wake' in f:
                    target_group = f['wake']
                
                keys = list(target_group.keys())
                real_timesteps = []
                for k in keys:
                    if k.startswith('timestep_'):
                        try:
                            real_timesteps.append(int(k.split('_')[1]))
                        except:
                            pass
                
                if real_timesteps:
                    real_max = max(real_timesteps)
                    print(f"HDF5 Scan: Found {len(real_timesteps)} steps. Max: {real_max} (Overriding config guess {max_step})")
                    max_step = real_max
                else:
                    print(f"HDF5 Scan: No timesteps found in {H5_PATH}.")
        except Exception as h5e:
            print(f"Warning: Failed to scan HDF5 file: {h5e}")

    # Turbine Params (for Diameter and dt)
    turb_config = config_data.get('turbine', {})
    model_name = turb_config.get('model', 'NTNU')
    wind_speed = turb_config.get('windSpeed', 10.0)
    tsr = turb_config.get('tsr', 6.0)
    
    # Load turbine_params.json
    turb_params_path = os.path.join(PROJECT_ROOT, 'data', model_name, 'turbine_params.json')
    with open(turb_params_path, 'r') as f:
        turb_data = json.load(f)
    r_tip = turb_data.get('rTip', 0.447)
    D = 2.0 * r_tip
    
    # Calculate dt
    if r_tip > 1e-6:
        omega = tsr * wind_speed / r_tip
        dt = (2.0 * np.pi / omega) / steps_per_rev
    else:
        dt = 0.01 # Fallback
        
    print(f"Loaded Turbine: {model_name}, D={D:.3f} m, dt={dt:.5f} s")
    
except Exception as e:
    print(f"Warning: Could not read config for defaults ({e}). using hardcoded fallback.")
    max_step = 1000
    D = 0.894 # Fallback NTNU default
    dt = 0.005 # Fallback

# Apply Logic
VORTICITY_TIMESTEP = args.t_inst if args.t_inst is not None else max_step

AVG_STEP = args.avg_step
AVG_END_T = args.t_end if args.t_end is not None else max_step
# Default start: 10 groups back (e.g. if step=10, 100 steps back)
AVG_START_T = args.t_start if args.t_start is not None else (AVG_END_T - 10 * AVG_STEP)
if AVG_START_T < 0: AVG_START_T = 0

# --- Configuration Mapping ---
# Plot Settings (None means auto-scale)
VORTICITY_VMIN = args.w_min
VORTICITY_VMAX = args.w_max
VELOCITY_VMIN = args.u_min
VELOCITY_VMAX = args.u_max

X_LIMITS = args.xlim
Y_LIMITS = args.ylim

# Analysis Parameters
# Prioritize absolute resolution if given, else calculate from PPD
if args.res_vort is not None:
    VORTICITY_RES = args.res_vort
else:
    VORTICITY_RES = D / args.ppd_vort
    print(f"Using Relative Resolution (Vorticity): {args.ppd_vort} PPD -> {VORTICITY_RES:.4f} m")

SLICE_MODE = args.slice_mode

if args.res_avg is not None:
    AVG_RES = args.res_avg
else:
    AVG_RES = D / args.ppd_avg
    print(f"Using Relative Resolution (Average): {args.ppd_avg} PPD -> {AVG_RES:.4f} m")

# Helper to run C++ tool
def run_postprocess(timestep, resolution, output_vtk):
    # Pass x_start/x_end factors to C++ tool (indices 7, 8)
    # We add a buffer to ensuring the grid covers the requested plot area
    # Buffer = 0.5D
    x_start_fac = X_LIMITS[0] - 0.5
    x_end_fac = X_LIMITS[1] + 0.5
    
    if not os.path.exists(EXECUTABLE):
        raise RuntimeError(f"Executable not found: {EXECUTABLE}")

    cmd = [
        EXECUTABLE,
        H5_PATH,
        output_vtk,
        CONFIG_PATH,
        str(timestep),
        str(resolution),
        SLICE_MODE,
        str(x_start_fac),
        str(x_end_fac)
    ]
    
    # Environment setup (copied from previous known working env)
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/apps/system/easybuild/software/libarchive/3.5.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/XZ/5.2.5-GCCcore-11.2.0/lib:/apps/system/easybuild/software/cURL/7.78.0-GCCcore-11.2.0/lib:/apps/system/easybuild/software/OpenSSL/1.1/lib:/apps/system/easybuild/software/bzip2/1.0.8-GCCcore-11.2.0/lib:/apps/system/easybuild/software/ncurses/6.2-GCCcore-11.2.0/lib:/apps/system/easybuild/software/HDF5/1.14.0-GCC-11.2.0-serial/lib:/apps/system/easybuild/software/Szip/2.1.1-GCCcore-11.2.0/lib:/apps/system/easybuild/software/binutils/2.37-GCCcore-11.2.0/lib:/apps/system/easybuild/software/zlib/1.2.11-GCCcore-11.2.0/lib:/apps/system/easybuild/software/GCCcore/11.2.0/lib64"
    
    # print(f"Processing t={timestep} (res={resolution}m)...")
    # Ensure C++ tool runs from PROJECT_ROOT so it can find data/ files relative to it
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=PROJECT_ROOT)
    if result.returncode != 0:
        print(f"Error at t={timestep}: {result.stderr}")
        raise RuntimeError("postprocess_grid failed")
    return output_vtk

# --- 1. Instantaneous Analysis (Vorticity & Velocity Deficit) ---
def analyze_instantaneous():
    tasks = args.tasks
    if 'all' in tasks:
        do_vor_z = True
        do_vor_mag = True
        do_inst_def = True
    else:
        do_vor_z = 'vor_z' in tasks
        do_vor_mag = 'vor_mag' in tasks
        do_inst_def = 'inst_def' in tasks
        
    if not (do_vor_z or do_vor_mag or do_inst_def):
        return

    print("\n--- Analyzing Instantaneous Fields ---")
    vtk_path = os.path.join(OUTPUT_DIR, f'inst_field_t{VORTICITY_TIMESTEP}.vtk')
    
    # Run postprocess once
    run_postprocess(VORTICITY_TIMESTEP, VORTICITY_RES, vtk_path)
    
    grid = pv.read(vtk_path)
    points = grid.points
    x = points[:, 0]
    y = points[:, 1]
    
    # Common Plot Helper
    def plot_scalar(scalar_data, filename, title, label, vmin=None, vmax=None, cmap='viridis'):
        plt.figure(figsize=(12, 4))
        
        if vmin is not None and vmax is not None:
            levels = np.linspace(vmin, vmax, 100)
        else:
            levels = 100
            
        plt.tricontourf(x/D, y/D, scalar_data, levels=levels, cmap=cmap, extend='both')
        
        cbar_ticks = None
        if vmin is not None and vmax is not None:
             cbar_ticks = np.linspace(vmin, vmax, 5)

        plt.colorbar(label=label, ticks=cbar_ticks)
        
        time_sec = VORTICITY_TIMESTEP * dt
        plt.title(f'{title}\nt={time_sec:.4f}s (step {VORTICITY_TIMESTEP})')
        plt.xlabel('x/D [-]')
        plt.ylabel('y/D [-]')
        plt.axis('scaled')
        
        if X_LIMITS is not None: plt.xlim(X_LIMITS)
        if Y_LIMITS is not None: plt.ylim(Y_LIMITS)
            
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)
        plt.close()
        print(f"Saved: {filename}")

    # 1. Vorticity (Z and Mag)
    if do_vor_z or do_vor_mag:
        # Compute vorticity vector directly
        grid = grid.compute_derivative(scalars="Velocity", vorticity=True)
        vorticity = grid['vorticity'] # (N, 3) vector
        
        if do_vor_z:
            omega_z = vorticity[:, 2] # Z component
            plot_scalar(omega_z, 'vorticity_z.png', 'Instantaneous Vertical Vorticity', 
                        'Vertical Vorticity $\omega_z$ [1/s]', VORTICITY_VMIN, VORTICITY_VMAX, 'RdBu_r')
            
        if do_vor_mag:
            omega_mag = np.linalg.norm(vorticity, axis=1) # Magnitude
            # Use positive range for magnitude. Default 0 to VMAX (or 300)
            mag_min = 0
            mag_max = VORTICITY_VMAX if VORTICITY_VMAX else 300
            plot_scalar(omega_mag, 'vorticity_mag.png', 'Instantaneous Vorticity Magnitude', 
                        'Vorticity Magnitude $|\omega|$ [1/s]', mag_min, mag_max, 'inferno')

    # 2. Velocity Deficit
    if do_inst_def:
        vel = grid['Velocity']
        ux = vel[:, 0]
        deficit = (wind_speed - ux) / wind_speed
        
        # Default range 0-0.8 or user specified
        def_min = VELOCITY_VMIN if VELOCITY_VMIN is not None else 0.0
        def_max = VELOCITY_VMAX if VELOCITY_VMAX is not None else 1.0
        plot_scalar(deficit, 'inst_deficit.png', 'Instantaneous Velocity Deficit', 
                    'Velocity Deficit [-]', def_min, def_max, 'viridis')

# --- 2. Time-averaged Streamwise Velocity ---
def analyze_time_average():
    tasks = args.tasks
    if not ('all' in tasks or 'avg_def' in tasks):
        return

    print("\n--- Analyzing Time-Averaged Velocity Deficit ---")
    timesteps = list(range(AVG_START_T, AVG_END_T + 1, AVG_STEP))
    
    sum_ux = None
    count = 0
    
    points_x = None
    points_y = None
    
    for t in timesteps:
        print(f"Processing t={t}...")
        vtk_path = os.path.join(OUTPUT_DIR, f'temp_vel_t{t}.vtk')
        
        # Always run to ensure fresh data (unless we want to cache)
        # Overwrite mode: always run
        run_postprocess(t, AVG_RES, vtk_path)
        
        grid = pv.read(vtk_path)
        vel = grid['Velocity'] # (N, 3)
        ux = vel[:, 0]
        
        if sum_ux is None:
            sum_ux = np.zeros_like(ux)
            points_x = grid.points[:, 0]
            points_y = grid.points[:, 1]
        
        sum_ux += ux
        count += 1
        
        # Clean up temp file
        os.remove(vtk_path)
        
    avg_ux = sum_ux / count
    
    # Calculate Velocity Deficit: (U_inf - U) / U_inf
    deficit = (wind_speed - avg_ux) / wind_speed
    
    # Plotting
    plt.figure(figsize=(12, 4))
    
    # Default range for deficit if not specified
    vmin = VELOCITY_VMIN if VELOCITY_VMIN is not None else 0.0
    vmax = VELOCITY_VMAX if VELOCITY_VMAX is not None else 1.0
    levels = np.linspace(vmin, vmax, 100)

    plt.tricontourf(points_x/D, points_y/D, deficit, levels=levels, cmap='viridis', extend='both')
    
    cbar_ticks = None
    if VELOCITY_VMIN is not None and VELOCITY_VMAX is not None:
        cbar_ticks = np.linspace(VELOCITY_VMIN, VELOCITY_VMAX, 5)
        
    plt.colorbar(label='Velocity Deficit $(U_{\infty} - U_x)/U_{\infty}$ [-]', ticks=cbar_ticks)
    
    t_start_sec = AVG_START_T * dt
    t_end_sec = AVG_END_T * dt
    plt.title(f'Time-Averaged Velocity Deficit\nt={t_start_sec:.2f}-{t_end_sec:.2f}s (step {AVG_START_T}-{AVG_END_T})')
    plt.xlabel('x/D [-]')
    plt.ylabel('y/D [-]')
    plt.axis('scaled')
    
    # View Control
    if X_LIMITS is not None:
        plt.xlim(X_LIMITS)
    if Y_LIMITS is not None:
        plt.ylim(Y_LIMITS)
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'mean_velocity_plot.png'), dpi=300)
    plt.close()
    print("Time-averaged velocity plot saved.")

if __name__ == "__main__":
    analyze_instantaneous()
    analyze_time_average()
    print("\nAnalysis Complete. Results in:", OUTPUT_DIR)
