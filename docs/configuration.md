# Configuration Reference (config.json)

FVW-CPP uses `config.json` to control simulation parameters, turbine models, and perturbation cases. This file is the single source of truth for your simulation settings.

## 1. Turbine Configuration

The `turbine` section defines the physical properties of the wind turbine and the inflow conditions.

```json
"turbine": {
    "model": "NREL_5MW",  // Turbine model name. Automatically loads geometry from data/NREL_5MW/
    "windSpeed": 11.4,    // Free stream wind speed (m/s)
    "nSegments": 18,      // Number of radial spanwise segments per blade
    "tsr": 7.0,           // Tip-Speed Ratio (lambda = omega * R / Uinf)
    "rho": 1.225          // Air density (kg/m^3), default is 1.225 if omitted
    
    // Optional Overrides:
    // If "model" is set, these are loaded from data/<model>/turbine_params.json.
    // You can override them here if needed:
    // "rTip": 63.0,      // Rotor radius (m)
    // "rHub": 1.5,       // Hub radius (m)
    // "nBlades": 3,      // Number of blades
    // "hubHeight": 90.0  // Hub height (m)
}
```

### Adding New Models
To add a new turbine model (e.g., `NTNU`):
1. Create `data/NTNU/turbine_params.json`.
2. Add blade geometry to `data/NTNU/blade_geometry.csv`.
3. Set `"model": "NTNU"` in your config.

## 2. Simulation Settings

The `simulation` section controls the time-stepping and numerical parameters of the Free Vortex Wake method.

```json
"simulation": {
    "dt": 0.05,               // Time step size (s)
    "totalTime": 10.0,        // Total simulation duration (s)
    "outputFrequency": 1,     // Write H5 output every N steps (1 = every step)
    "cutoffParam": 0.1,       // Vortex core cut-off parameter (regularization)
    "coreType": "ChordBasedCore", // Core model: "VanGarrel" or "ChordBasedCore"
    "vortexModel": "Constant"     // Vortex decay: "Constant" or "GammaDecay"
}
```

## 3. Perturbation Cases

You can define multiple simulation cases in the `perturbations` array. The main program runs them sequentially.

```json
"perturbations": [
    {
        "name": "case_baseline",   // Output folder name (results/case_baseline/)
        "type": "None",            // Perturbation type
        "amplitude": 0.0,
        "freqFactor": 0.0
    },
    {
        "name": "case_dynamic_pitch",
        "type": "CollectivePitch", 
        "amplitude": 1.0,          // Amplitude in degrees
        "freqFactor": 1.5          // Frequency relative to rotor speed (f = 1.5 * f_rotor)
    }
]
```

### Supported Perturbation Types
| Type | Description | Parameters |
|------|-------------|------------|
| `None` | Baseline constant operation | None |
| `CollectivePitch` | Dynamic sinusoidal collective pitch | `amplitude` (deg), `freqFactor` |
| `AsymmetricStaticPitch` | Static offset applied to one blade | `amplitude` (deg) |
