
# Research Plan: Tip Vortex Instability Mechanisms

## 1. Objective
Investigate the onset and evolution of tip vortex instabilities (specifically **Crow Instability** and **Pairing/Leapfrogging**) using the FVW-CPP solver.
The goal is to quantify how parameters like **Tip Speed Ratio (TSR)**, **Vortex Core Size (Cutoff)**, and **Perturbation Amplitude** affect the instability growth rate and breakdown location.

## 2. Methodology
The study utilizes the Inviscid Free Vortex Wake (FVW) method. Since the model is naturally stable, **external perturbations** must be introduced to trigger physically meaningful instabilities.

### Key Variables
*   **Cutoff Parameter (`cutoffParam`)**: Controls the vortex core radius.
    *   *Range*: 0.05 (Fine) to 0.25 (Coarse).
    *   *Hypothesis*: Smaller cutoff $\to$ Faster short-wave growth (Pairing).
*   **Perturbation**:
    *   *Type*: `AsymmetricStaticPitch` (Geometry) or `CollectivePitch` (Transient).
    *   *Amplitude*: $0.5^\circ \sim 2.0^\circ$.
*   **TSR**: Tip Speed Ratio (Affects helix pitch and mutual inductance).

## 3. Workflow & Tools
All analysis scripts are located in `tools/python/`.

### Step 1: Simulation Setup
Edit `config_ntnu.json`:
```json
"perturbations": [
    {
        "type": "AsymmetricStaticPitch", 
        "amplitude": 1.0, 
        "freqFactor": 0.0
    }
],
"simulation": {
    "cutoffParam": 0.1
}
```

### Step 2: Execution
Run the simulation on the cluster:
```bash
sbatch submit_job.slurm
```
(Ensure `outputFrequency` is small enough, e.g., 10, to capture dynamics).

### Step 3: Post-Processing & Analysis

#### A. Instability Analysis (Primary Tool)
Quantify radial wobble and pairing.
*   **Script**: `tools/python/analyze_instability.py`
*   **Usage**:
    ```bash
    python tools/python/analyze_instability.py
    ```
*   **Outputs** (`results/NTNU_Baseline/post_processing/instability_analysis/`):
    *   `radial_instability.png`: Visualizes $r/x$ trajectory. Look for oscillating growth (Crow) or straight lines (Stable).
    *   `pairing_distance.png`: Distance between adjacent vortices. Look for sudden drops (Leapfrogging).

#### B. Validation with ALM-LES
Compare mean velocity profiles with high-fidelity LES data.
*   **Script**: `tools/python/compare_fvw_alm.py`
*   **Usage**:
    ```bash
    python tools/python/compare_fvw_alm.py
    ```
*   **Requires**: `results/NTNU_Baseline/probe_output.csv` (Run simulation with `computeProbes: true`).
*   **Outputs**: `results/NTNU_Baseline/post_processing/comparison_plots/`

#### C. Wake Metrics
Compute deficit and expansion trends.
*   **Script**: `tools/python/analyze_wake_metrics.py` (Optional)

## 4. File Cleanup Strategy
To save space and keep the workspace clean:
*   **Keep**: `wake.h5` (Core data), `config_*.json` (Record of settings), `analysis_*` (Plots).
*   **Delete**: 
    *   `probe_output_partial.csv` (Intermediate debug file).
    *   `slurm-*.out` (Old logs, unless containing errors).
    *   `debug_*.txt` (If any).

## 5. Next Steps Experiment Matrix

| Experiment ID | Cutoff | Perturbation | Objective |
| :--- | :--- | :--- | :--- |
| **Exp-01 (Baseline)** | 0.1 | None | Confirm absolute stability (Straight lines). |
| **Exp-02 (Low Noise)**| 0.25 | None | Verify removal of numerical noise/ejections. |
| **Exp-03 (Forced)**   | 0.25 | 1.0 deg Pitch | Trigger & Measure Crow Instability Growth Rate. |
| **Exp-04 (Fine Core)**| 0.05 | 1.0 deg Pitch | Investigate Short-wave (Pairing) sensitivity. |

