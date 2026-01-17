import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import h5py

def main():
    parser = argparse.ArgumentParser(description="FVW Wake ZX Projection (Side View)")
    parser.add_argument("file", help="Path to wake.h5")
    parser.add_argument("--step", type=int, default=-1, help="Timestep to analyze (default: last step)")
    parser.add_argument("--blade", type=int, default=0, help="Blade index to analyze (default: 0)")
    args = parser.parse_args()

    h5_path = args.file
    if not os.path.exists(h5_path):
        print(f"Error: File {h5_path} not found.")
        return

    print(f"==================================================")
    print(f"   FVW Wake ZX Projection")
    print(f"   File: {h5_path}")
    print(f"==================================================")

    try:
        f = h5py.File(h5_path, "r")
        
        # Get Timestep
        timesteps = []
        if 'wake' in f:
             for k in f['wake'].keys():
                if k.startswith("timestep_"):
                    timesteps.append(int(k.split("_")[-1]))
             timesteps.sort()
        
        if not timesteps:
            print("Error: No timesteps found.")
            return

        t_idx = args.step
        if t_idx == -1:
            t_idx = timesteps[-1]
        
        print(f"Analyzing Timestep: {t_idx}")
        
        # Read Blade Data
        grp = f[f'wake/timestep_{t_idx}/blade_{args.blade}']
        nodes = np.array(grp['nodes'][:]) # (N, 6)
        lines = np.array(grp['lines'][:]) # (M, 4) -> start, end, gamma, type
        
        # Filter Logic:
        # User wants "Outer Filament Data" (Tip Vortex)
        # And projected on ZX plane (Side View), filtering for Y > 0 (to avoid clutter from the other side if rotating?)
        # But for a single blade analysis, we might want to see everything or just the tip.
        
        # Let's perform the Robust Tip Vortex Trace first to isolate the filament
        # 1. Build Adjacency
        adj_map = {} 
        for i in range(len(lines)):
            if lines[i, 3] == 2.0: # Trailing
                end_idx = int(lines[i, 1])
                if end_idx not in adj_map:
                    adj_map[end_idx] = []
                adj_map[end_idx].append(i) 

        # 2. Identify Tip Node (Max Radius)
        potential_tips = list(adj_map.keys())
        max_r = -1.0
        tip_node_idx = -1
        for idx in potential_tips:
            pos = nodes[idx, :3]
            r = np.sqrt(pos[1]**2 + pos[2]**2)
            if r > max_r:
                max_r = r
                tip_node_idx = idx
        
        if tip_node_idx == -1:
             print("Error: Tip node not found.")
             return

        # 3. Trace Filament
        filament_nodes = []
        filament_gammas = []
        
        curr_node = tip_node_idx
        # filament_nodes.append(nodes[curr_node, :3]) # Start point
        # filament_gammas.append(np.nan)
        
        while curr_node in adj_map:
            line_idx = adj_map[curr_node][0]
            line_data = lines[line_idx]
            start_node = int(line_data[0])
            gamma = line_data[2]
            
            p = nodes[start_node, :3]
            filament_nodes.append(p)
            filament_gammas.append(gamma)
            
            curr_node = start_node
            if len(filament_nodes) > 10000: break
            
        filament_nodes = np.array(filament_nodes)
        filament_gammas = np.array(filament_gammas)
        
        if len(filament_nodes) == 0:
            print("No filament traced.")
            return

        # 4. Plot ZX Projection (Side View)
        # Filter for Y > 0 as requested by user logic (though for single filament trace, it might not matter unless it spirals)
        # Actually, for a single tip vortex, filtering y>0 slices the helix. 
        # If we plot distinct points where y>0, we see the "top" of the loops.
        
        # Let's replicate the user's "Slice" effect
        # The user code: mask = (pos[:, 0] >= 0.0) & (pos[:, 1] > 0.0)
        # This keeps only the "Upper" (or one side) half of the helical turns.
        
        pos = filament_nodes
        gamma = filament_gammas
        
        mask = (pos[:, 0] >= 0.0) & (pos[:, 1] > 0.0)
        pos_filtered = pos[mask]
        gamma_filtered = gamma[mask]
        
        fig, ax = plt.subplots(figsize=(12, 5))
        
        sc = ax.scatter(pos_filtered[:, 0], pos_filtered[:, 2], c=gamma_filtered, cmap='viridis', s=20)
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label('Circulation $\Gamma$')
        
        ax.set_xlabel('Downstream Distance X [m]')
        ax.set_ylabel('Vertical Position Z [m]')
        ax.set_title(f"Tip Vortex Side View (ZX Projection, Y>0) - t={t_idx}")
        ax.grid(True, alpha=0.5, linestyle='--')
        ax.set_aspect('equal')
        
        plt.tight_layout()
        
        out_dir = os.path.dirname(h5_path)
        out_png = os.path.join(out_dir, f"zx_projection_t{t_idx}.png")
        plt.savefig(out_png, dpi=150)
        print(f"\n[Success] Plot saved to: {os.path.abspath(out_png)}")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
