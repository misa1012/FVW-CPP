import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import h5py

def main():
    parser = argparse.ArgumentParser(description="FVW Wake Visualization Tool")
    parser.add_argument("file", help="Path to wake.h5")
    parser.add_argument("--step", type=int, default=-1, help="Timestep to analyze (default: last step)")
    parser.add_argument("--save", action="store_true", help="Save as PNG instead of showing")
    args = parser.parse_args()

    h5_path = args.file
    if not os.path.exists(h5_path):
        print(f"Error: File {h5_path} not found.")
        return

    print(f"==================================================")
    print(f"   FVW Wake Visualization")
    print(f"   File: {h5_path}")
    print(f"==================================================")

    try:
        f = h5py.File(h5_path, "r")
        timesteps = []
        if 'wake' in f:
             # Find available timesteps
             for k in f['wake'].keys():
                if k.startswith("timestep_"):
                    try:
                        timesteps.append(int(k.split("_")[-1]))
                    except ValueError:
                        pass
             timesteps.sort()
        else:
            print("Error: No 'wake' group in HDF5 file.")
            return
            
        if not timesteps:
            print("Error: No timesteps found.")
            return

        # Select Timestep
        t_idx = args.step
        if t_idx == -1:
            t_idx = timesteps[-1]
        
        if t_idx not in timesteps:
            print(f"Error: Timestep {t_idx} not found.")
            return

        print(f"Visualizing Timestep: {t_idx}")
        
        # Read Geometry (Blades)
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        grp = f[f'wake/timestep_{t_idx}']
        
        colors = ['b', 'r', 'g']
        
        # Iterate over blades
        n_blades = 0
        for key in grp.keys():
            if key.startswith("blade_"):
                b_idx = int(key.split("_")[-1])
                nodes = np.array(grp[f'{key}/nodes'][:]) # (N, 6) -> x,y,z, u,v,w
                
                # Plot X, Y, Z (Columns 0, 1, 2)
                # Scatter points
                # ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], s=5, c=colors[b_idx % len(colors)])
                
                # We need to reconstruct lines. 
                # Assuming nodes are ordered: Spanwise first, then Trailing? 
                # Actually, nodes are unstructured in a simple list.
                # But 'lines' dataset tells connectivity.
                
                # Quick hack: Just plot points for now to verify geometry
                ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], s=2, c=colors[b_idx % len(colors)], label=f'Blade {b_idx}')
                n_blades += 1

        # Set labels and limits
        ax.set_xlabel('X (Axial)')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f"Wake Geometry at t={t_idx}")
        
        # Try to set equal aspect ratio
        try:
             # Create cubic bounding box to force equal aspect ratio
             x_limits = ax.get_xlim3d()
             y_limits = ax.get_ylim3d()
             z_limits = ax.get_zlim3d()

             x_range = abs(x_limits[1] - x_limits[0])
             x_middle = np.mean(x_limits)
             y_range = abs(y_limits[1] - y_limits[0])
             y_middle = np.mean(y_limits)
             z_range = abs(z_limits[1] - z_limits[0])
             z_middle = np.mean(z_limits)

             # The plot bounding box is a sphere in the sense of the infinite
             # teem, but using box for this.
             plot_radius = 0.5*max([x_range, y_range, z_range])

             ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
             ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
             ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
        except Exception:
             pass

        ax.legend()
        
        if args.save:
             out_dir = os.path.dirname(h5_path)
             out_png = os.path.join(out_dir, f"wake_viz_t{t_idx}.png")
             plt.savefig(out_png, dpi=150)
             print(f"[Success] Saved to {out_png}")
        else:
             print("Saving to PNG by default in non-interactive mode")
             out_dir = os.path.dirname(h5_path)
             out_png = os.path.join(out_dir, f"wake_viz_t{t_idx}.png")
             plt.savefig(out_png, dpi=150)
             print(f"[Success] Saved to {out_png}")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
