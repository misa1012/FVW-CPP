import h5py
import numpy as np
import os

class WakeReader:
    """
    Handles reading of FVW-CPP HDF5 output files.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
    def list_timesteps(self):
        """Returns sorted list of available timesteps."""
        with h5py.File(self.filepath, "r") as f:
            if "/liftingline" not in f:
                return []
            ts = []
            for k in f["/liftingline"].keys():
                if k.startswith("timestep_"):
                    try:
                        ts.append(int(k.split("_")[-1]))
                    except ValueError:
                        pass
            return sorted(ts)

    def read_config(self):
        """
        Reads simulation configuration and geometry.
        Returns a dictionary containing simulation parameters.
        """
        config = {}
        with h5py.File(self.filepath, "r") as f:
            # Simulation params
            cfg_sim = f.get("/config/simulation")
            if cfg_sim:
                config['rho'] = float(cfg_sim['rho'][0])
                config['Uinf'] = float(cfg_sim['wind_speed'][0])
                config['omega'] = float(cfg_sim['omega'][0])
                config['dt'] = float(cfg_sim['dt'][0])
                config['n_blades'] = int(cfg_sim['n_blades'][0])
                # Some files might call it r_tip or R
                if 'r_tip' in cfg_sim:
                    config['r_tip'] = float(cfg_sim['r_tip'][0])
                else:
                    config['r_tip'] = 63.0 # Fallback or infer

            # Geometry
            cfg_geom = f.get("/config/geometry")
            if cfg_geom:
                config['r_center'] = np.array(cfg_geom['r_shed'][:])
                config['r_bounds'] = np.array(cfg_geom['r_trail'][:])
                config['chord'] = np.array(cfg_geom['chord'][:])
                config['dr'] = np.diff(config['r_bounds'])
                
                if 'twist' in cfg_geom:
                    config['twist'] = np.array(cfg_geom['twist'][:])
                else:
                    # Fallback for NREL 5MW if not saved (legacy support)
                    config['twist'] = self._get_fallback_twist(config['r_bounds'], config['r_center'])
        
        return config

    def _get_fallback_twist(self, r_bounds, r_center):
        """Fallback twist distribution for NREL 5MW."""
        original_r = np.array([1.5, 2.86670, 5.6, 8.33330, 11.75, 15.85, 19.95,
                               24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65,
                               52.75, 56.1667, 58.90, 61.6333, 63.0])
        original_twist = -1 * np.array([13.308, 13.308, 13.308, 13.308, 13.308, 11.48,
                                        10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125,
                                        2.319, 1.526, 0.863, 0.370, 0.106, 0.106/2])
        twist_trailing = np.interp(r_bounds, original_r, original_twist)
        return np.interp(r_center, r_bounds, twist_trailing)

    def read_timestep_data(self, timestep):
        """
        Reads performance and velocity data for a specific timestep.
        Returns:
            list of dicts, one per blade, with keys: 'cl', 'cd', 'aoa', 'rel_vel'
        """
        data = []
        with h5py.File(self.filepath, "r") as f:
            base_path = f"/liftingline/timestep_{timestep}"
            if base_path not in f:
                return None
            
            # Assuming n_blades is consistent, check keys
            blade_keys = [k for k in f[base_path].keys() if k.startswith("blade_")]
            n_blades = len(blade_keys)
            
            for b in range(n_blades):
                blade_path = f"{base_path}/blade_{b}"
                perf = np.array(f[f"{blade_path}/perf"][:]) # [nSeg, 3] -> cl, cd, aoa
                vel = np.array(f[f"{blade_path}/relative_velocity_bcs"][:]) # [nSeg, 3]
                
                data.append({
                    'cl': perf[:, 0],
                    'cd': perf[:, 1],
                    'aoa': perf[:, 2],
                    'rel_vel': vel
                })
        return data
