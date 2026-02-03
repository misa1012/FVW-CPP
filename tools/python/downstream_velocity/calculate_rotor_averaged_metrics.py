
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Rotor-Averaged Metrics (REWS, TI) from ALM Slices.")
    parser.add_argument("--alm_dir", required=True, help="Path to ALM-LES case directory")
    parser.add_argument("--time", default="1.872", help="Time step to analyze (default: 1.872)")
    parser.add_argument("--u_inf", type=float, default=10.0, help="Freesteam velocity (default: 10.0)")
    parser.add_argument("--diameter", type=float, default=0.894, help="Rotor Diameter D (default: 0.894)")
    parser.add_argument("--hub_height", type=float, default=0.8, help="Hub Height (default: 0.8)")
    parser.add_argument("--output_dir", default="results/alm_metrics", help="Output directory for plots and data")
    return parser.parse_args()

def calculate_metrics_for_slice(vtk_path, fvw_D, hub_height, U_inf):
    """
    Calculates REWS and TI for a single slice file.
    Assumes fields: UMean, UPrime2Mean (or kMean)
    """
    try:
        grid = pv.read(vtk_path)
    except Exception as e:
        print(f"Error reading {vtk_path}: {e}")
        return None

    # Points
    pts = grid.points
    y = pts[:, 1]
    z = pts[:, 2]
    
    # Calculate distance from hub center
    r_dist = np.sqrt(y**2 + (z - hub_height)**2)
    R = fvw_D / 2.0
    
    # Filter for rotor disk
    mask = r_dist <= R
    
    if not np.any(mask):
        return None
    
    # Extract Data within Mask
    # Check available fields
    field_names = grid.point_data.keys()
    
    # 1. Axial Velocity (UMean or U)
    # Prefer UMean if available, else U
    u_field_name = "UMean" if "UMean" in field_names else "U"
    if u_field_name not in field_names:
        print(f"Warning: Neither UMean nor U found in {vtk_path}")
        return None
        
    U_vec = grid.point_data[u_field_name][mask]
    Ux = U_vec[:, 0] # Axial component
    
    # 2. Turbulence Intensity (TI)
    # TI = sqrt(2/3 * k) / U_inf (or local U, standard usually U_hub or U_inf)
    # Need 'kMean' or 'UPrime2Mean'
    k = None
    if "turbulenceProperties:kMean" in field_names:
        k = grid.point_data["turbulenceProperties:kMean"][mask]
        print(f"Using turbulenceProperties:kMean for TI")
    elif "kMean" in field_names:
        k = grid.point_data["kMean"][mask]
        print(f"Using kMean for TI")
    elif "UPrime2Mean" in field_names:
        # Reynolds stress tensor: xx xy xz yy yz zz
        # k = 0.5 * (uu + vv + ww)
        # In OpenFOAM, symmTensor usually stored as: xx xy xz yy yz zz
        # Components: 0:xx, 1:xy, 2:xz, 3:yy, 4:yz, 5:zz
        up2 = grid.point_data["UPrime2Mean"][mask]
        # Check shape, usually (N, 6) for symmTensor
        if up2.shape[1] == 6:
            uu = up2[:, 0]
            vv = up2[:, 3]
            ww = up2[:, 5]
            k = 0.5 * (uu + vv + ww)
            print(f"Calculated k from UPrime2Mean")
        else:
             print(f"UPrime2Mean has unexpected shape {up2.shape}")

    # --- Calculations ---
    
    # A. REWS
    # REWS = cuberoot( 1/A * integral(U^3 dA) )
    # Since points are unstructured, simple mean might be biased if density varies.
    # However, 'surfaces' usually outputs fairly uniform triangulation or based on mesh.
    # For now, we assume simple area-weighted mean is approximated by arithmetic mean if uniform.
    # OpenFOAM 'cuttingPlane' usually generates triangles. PyVista can compute cell areas.
    
    # Better approach: Use cell areas for weighting if possible.
    # grid is PolyData.
    # Let's try simple mean first (often sufficient for 'uniform' source mesh)
    
    u3_mean = np.mean(Ux**3)
    rews = np.cbrt(u3_mean)
    rews_norm = rews / U_inf
    
    # B. Rotor-Averaged TI
    # TI_local = sqrt(2/3 * k) / U_local (or U_inf)
    # Usually 'Turbulence Intensity' implies normalization by U_hub or U_inf.
    # The reference plot shows TI ~ 0.1-0.4.
    # Formula: TI_avg = 1/A * integral( TI_local dA ) -> No
    # Usually it is: TI_avg = sqrt( 2/3 * k_avg ) / U_inf   OR   avg( sqrt(2/3*k)/U_local )
    # Standard IEC: TI = sigma_u / U_hub.
    # For CFD (k-based): TI = sqrt(2/3 * k) / U_ref.
    # Let's compute average k first.
    
    ti_avg = np.nan
    if k is not None:
        k_mean_val = np.mean(k)
        # Using U_inf as reference velocity is standard for "TI vs Distance"
        # If using local velocity, singularity at low U.
        ti_avg = np.sqrt(2.0/3.0 * k_mean_val) / U_inf
    
    return {
        "REWS_norm": rews_norm,
        "TI_avg": ti_avg
    }

def main():
    args = parse_args()
    
    # Directory with slices
    # We expect 'surfaces_yz_highres'
    vtk_dir = os.path.join(args.alm_dir, "postProcessing", "surfaces_yz_highres", args.time)
    
    if not os.path.exists(vtk_dir):
        print(f"Wait! Directory not found: {vtk_dir}")
        print("Did you run postProcess with the highres dictionary?")
        return # Don't exit error, just validation
        
    os.makedirs(args.output_dir, exist_ok=True)
    
    import glob
    files = glob.glob(os.path.join(vtk_dir, "yz_*D.vtp")) + glob.glob(os.path.join(vtk_dir, "yz_*D.vtk"))
    
    results = []
    
    print(f"Processing {len(files)} slices...")
    
    for f in files:
        fname = os.path.basename(f)
        # Parse distance from filename: yz_neg2_0D.vtp or yz_3_5D.vtp
        base = fname.split('.')[0] # yz_...
        
        # Remove 'yz_' prefix
        dist_str = base.replace('yz_', '').replace('D', '')
        # Handle 'neg' -> '-' and '_' -> '.'
        dist_str = dist_str.replace('neg', '-').replace('_', '.')
        
        try:
            dist = float(dist_str)
        except:
            print(f"Skipping {fname}, cannot parse distance.")
            continue
            
        metrics = calculate_metrics_for_slice(f, args.diameter, args.hub_height, args.u_inf)
        if metrics:
            metrics["x/D"] = dist
            results.append(metrics)
            
    # Sort by distance
    results.sort(key=lambda x: x["x/D"])
    
    df = pd.DataFrame(results)
    print(df)
    
    # --- Plotting ---
    # 1. REWS
    plt.figure(figsize=(10, 5))
    plt.plot(df["x/D"], df["REWS_norm"], 'o-', label="ALM-LES")
    plt.xlabel("x/D")
    plt.ylabel("REWS / U_inf")
    plt.title("Rotor Equivalent Wind Speed (REWS)")
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(args.output_dir, "rews_vs_distance.png"), dpi=300)
    plt.close()
    
    # 2. Rotor Averaged TI
    plt.figure(figsize=(10, 5))
    plt.plot(df["x/D"], df["TI_avg"], 's-', color='orange', label="ALM-LES")
    plt.xlabel("x/D")
    plt.ylabel("Rotor Averaged TI")
    plt.title("Turbulence Intensity (Streamwise)")
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(args.output_dir, "ti_vs_distance.png"), dpi=300)
    plt.close()
    
    print(f"Saved plots to {args.output_dir}")

if __name__ == "__main__":
    main()
