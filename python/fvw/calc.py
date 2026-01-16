import numpy as np
from .io import WakeReader

class PerformanceCalculator:
    """
    Calculates aerodynamic performance metrics (Thrust, Power, Cp, Ct) 
    from simulation results.
    """
    def __init__(self, reader: WakeReader):
        self.reader = reader
        self.config = reader.read_config()
        self.timesteps = reader.list_timesteps()

    def compute_time_series(self, start_rev=0.0, end_time=None):
        """
        Computes integrated loads over time.
        
        Args:
            start_rev (float): Start calculation after these many rotor revolutions (to skip transient).
            end_time (float): Optional end time in seconds.
            
        Returns:
            dict: time, rev, CP, CT, Power, Thrust arrays.
        """
        cfg = self.config
        rho, U, omega, dt = cfg['rho'], cfg['Uinf'], cfg['omega'], cfg['dt']
        r_tip = cfg['r_tip']
        n_blades = cfg['n_blades']
        
        # Geometry
        r = cfg['r_center']
        dr = cfg['dr']
        chord = cfg['chord']
        twist = cfg['twist']
        
        # Reference Area
        A = np.pi * r_tip**2
        
        # Determine timestep range
        ts = self.timesteps
        if not ts:
            return None
            
        T_rev = 2 * np.pi / omega
        start_step = int(np.ceil(start_rev * T_rev / dt))
        
        selected_ts = [t for t in ts if t >= start_step]
        if end_time is not None:
            max_step = int(np.floor(end_time / dt))
            selected_ts = [t for t in selected_ts if t <= max_step]
            
        if not selected_ts:
            return None

        # Storage
        results = {
            'time': [],
            'CP': [],
            'CT': [],
            'Power': [],
            'Thrust': []
        }
        
        print(f"Processing {len(selected_ts)} timesteps...")
        
        # Iteration
        for t in selected_ts:
            blade_data = self.reader.read_timestep_data(t)
            
            torque_sum = 0.0
            thrust_sum = 0.0
            
            for b in range(n_blades):
                bd = blade_data[b]
                cl, cd, aoa = bd['cl'], bd['cd'], bd['aoa']
                rel_vel = bd['rel_vel']
                
                # Magnitude squared of relative velocity
                v2 = np.sum(rel_vel**2, axis=1)
                
                # Inflow angle phi = AoA - Twist (all in degrees in file? No, usually AoA deg, twist deg)
                # Note: The notebook example assumed: phi = radians(aoa - twist)
                # But typically AoA = phi - twist => phi = AoA + twist ???
                # Let's check the notebook logic again:
                # "phi=np.radians(aoa-twist)" 
                # This implies AoA output from code includes twist? Or is geometric AoA?
                # Actually commonly: alpha = phi - theta (twist). So phi = alpha + theta.
                # The notebook has: phi = aoa - twist. This is weird unless 'twist' is defined negative?
                # In notebook: "original_twist = -1 * np.array([...])". So twist is negative.
                # So phi = aoa - (-|twist|) = aoa + |twist|. Correct.
                
                phi = np.radians(aoa - twist)
                
                # Forces per unit length
                # L = 0.5 * rho * v^2 * c * Cl
                L = 0.5 * rho * v2 * chord * cl
                D = 0.5 * rho * v2 * chord * cd
                
                # Projetion
                # dT = (L cos phi + D sin phi) * dr
                # dQ = (L sin phi - D cos phi) * r * dr
                
                dT = (L * np.cos(phi) + D * np.sin(phi)) * dr
                dQ = (L * np.sin(phi) - D * np.cos(phi)) * r * dr
                
                thrust_sum += np.sum(dT)
                torque_sum += np.sum(dQ)
                
            power = torque_sum * omega
            
            # Coefficients
            # CT = T / (0.5 rho U^2 A)
            # CP = P / (0.5 rho U^3 A)
            
            CT = thrust_sum / (0.5 * rho * U**2 * A)
            CP = power / (0.5 * rho * U**3 * A)
            
            current_time = t * dt
            
            results['time'].append(current_time)
            results['CP'].append(CP)
            results['CT'].append(CT)
            results['Power'].append(power)
            results['Thrust'].append(thrust_sum)
            
        # Convert to numpy arrays
        for k in results:
            results[k] = np.array(results[k])
            
        results['rev'] = results['time'] / T_rev
        results['omega'] = omega
        
        return results
