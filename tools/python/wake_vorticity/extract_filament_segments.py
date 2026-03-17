#!/usr/bin/env python3
"""
Extract filament segments from FVW wake.h5 for later circulation-to-grid projection.

Outputs:
- compressed NPZ with segment geometry and metadata
- plain-text summary for quick inspection

The extracted segment representation is:
    center, tangent, length, gamma, line_type
where the circulation vector content of one segment is
    gamma * tangent * length
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


TYPE_LABELS = {
    0: "bound",
    1: "shed",
    2: "trailing",
}


def find_timesteps(h5: h5py.File) -> list[int]:
    if "wake" not in h5:
        raise RuntimeError("Missing /wake group in wake.h5")
    return sorted(
        int(k.split("_")[-1])
        for k in h5["wake"].keys()
        if k.startswith("timestep_")
    )


def choose_timestep(h5: h5py.File, timestep: str) -> int:
    timesteps = find_timesteps(h5)
    if not timesteps:
        raise RuntimeError("No /wake/timestep_* groups found")
    if timestep == "last":
        return timesteps[-1]
    t = int(timestep)
    if t not in timesteps:
        raise RuntimeError(f"Timestep {t} not found. Available range: {timesteps[0]}..{timesteps[-1]}")
    return t


def load_segments(h5_path: Path, timestep: int) -> dict[str, np.ndarray]:
    blade_ids = []
    start_xyz = []
    end_xyz = []
    center_xyz = []
    tangent_xyz = []
    lengths = []
    gammas = []
    line_types = []
    start_idx = []
    end_idx = []

    with h5py.File(h5_path, "r") as h5:
        grp = h5[f"/wake/timestep_{timestep}"]
        blade_names = sorted(k for k in grp.keys() if k.startswith("blade_"))
        if not blade_names:
            raise RuntimeError(f"No blade_* groups found in /wake/timestep_{timestep}")

        for blade_name in blade_names:
            blade_id = int(blade_name.split("_")[-1])
            blade_grp = grp[blade_name]
            nodes = np.asarray(blade_grp["nodes"][:, :3], dtype=float)
            lines = np.asarray(blade_grp["lines"][:], dtype=float)
            if lines.size == 0:
                continue

            i0 = lines[:, 0].astype(int)
            i1 = lines[:, 1].astype(int)
            gamma = lines[:, 2].astype(float)
            line_type = lines[:, 3].astype(int)

            p0 = nodes[i0]
            p1 = nodes[i1]
            d = p1 - p0
            ds = np.linalg.norm(d, axis=1)
            good = ds > 0.0
            if not np.all(good):
                i0 = i0[good]
                i1 = i1[good]
                gamma = gamma[good]
                line_type = line_type[good]
                p0 = p0[good]
                p1 = p1[good]
                d = d[good]
                ds = ds[good]

            t_hat = d / ds[:, None]
            c = 0.5 * (p0 + p1)

            n = len(ds)
            blade_ids.append(np.full(n, blade_id, dtype=int))
            start_xyz.append(p0)
            end_xyz.append(p1)
            center_xyz.append(c)
            tangent_xyz.append(t_hat)
            lengths.append(ds)
            gammas.append(gamma)
            line_types.append(line_type)
            start_idx.append(i0)
            end_idx.append(i1)

    if not lengths:
        raise RuntimeError(f"No valid segments found in {h5_path} at timestep {timestep}")

    return {
        "blade_id": np.concatenate(blade_ids),
        "start_xyz": np.vstack(start_xyz),
        "end_xyz": np.vstack(end_xyz),
        "center_xyz": np.vstack(center_xyz),
        "tangent_xyz": np.vstack(tangent_xyz),
        "length": np.concatenate(lengths),
        "gamma": np.concatenate(gammas),
        "line_type": np.concatenate(line_types),
        "start_idx": np.concatenate(start_idx),
        "end_idx": np.concatenate(end_idx),
    }


def write_summary(out_txt: Path, h5_path: Path, timestep: int, data: dict[str, np.ndarray]) -> None:
    blade_id = data["blade_id"]
    gamma = data["gamma"]
    length = data["length"]
    line_type = data["line_type"]
    tangent = data["tangent_xyz"]

    circulation_vector = gamma[:, None] * tangent * length[:, None]
    total_vec = circulation_vector.sum(axis=0)

    lines = []
    lines.append("Filament segment extraction summary")
    lines.append(f"h5_path: {h5_path}")
    lines.append(f"timestep: {timestep}")
    lines.append("")
    lines.append(f"n_segments: {len(length)}")
    lines.append(f"n_blades: {len(np.unique(blade_id))}")
    lines.append(f"total_length: {length.sum():.6f} m")
    lines.append(f"gamma_min: {gamma.min():.6e}")
    lines.append(f"gamma_max: {gamma.max():.6e}")
    lines.append(f"mean_abs_gamma: {np.mean(np.abs(gamma)):.6e}")
    lines.append(f"max_length: {length.max():.6e} m")
    lines.append(f"min_length: {length.min():.6e} m")
    lines.append("")
    lines.append("Integrated circulation-vector content sum(gamma * t_hat * ds):")
    lines.append(f"  x: {total_vec[0]:.6e}")
    lines.append(f"  y: {total_vec[1]:.6e}")
    lines.append(f"  z: {total_vec[2]:.6e}")
    lines.append("")
    lines.append("By line type:")
    for type_id in sorted(np.unique(line_type)):
        m = line_type == type_id
        label = TYPE_LABELS.get(int(type_id), f"unknown_{int(type_id)}")
        vec = circulation_vector[m].sum(axis=0)
        lines.append(
            f"  {label}: n={m.sum()}, length_sum={length[m].sum():.6f} m, "
            f"gamma_range=[{gamma[m].min():.6e}, {gamma[m].max():.6e}], "
            f"sum(gamma*t*ds)=({vec[0]:.6e}, {vec[1]:.6e}, {vec[2]:.6e})"
        )
    lines.append("")
    lines.append("By blade:")
    for b in sorted(np.unique(blade_id)):
        m = blade_id == b
        vec = circulation_vector[m].sum(axis=0)
        lines.append(
            f"  blade_{int(b)}: n={m.sum()}, length_sum={length[m].sum():.6f} m, "
            f"mean_abs_gamma={np.mean(np.abs(gamma[m])):.6e}, "
            f"sum(gamma*t*ds)=({vec[0]:.6e}, {vec[1]:.6e}, {vec[2]:.6e})"
        )

    out_txt.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract filament segments from FVW wake.h5")
    parser.add_argument("wake_h5", help="Path to wake.h5")
    parser.add_argument("--timestep", default="last", help="Wake timestep to extract, or 'last'")
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory. Default: <case>/post_processing/filament_segments",
    )
    args = parser.parse_args()

    h5_path = Path(args.wake_h5).resolve()
    if not h5_path.exists():
        raise FileNotFoundError(f"wake.h5 not found: {h5_path}")

    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else h5_path.parent / "post_processing" / "filament_segments"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as h5:
        timestep = choose_timestep(h5, args.timestep)

    data = load_segments(h5_path, timestep)

    stem = f"timestep_{timestep}"
    out_npz = out_dir / f"segments_{stem}.npz"
    out_txt = out_dir / f"segments_summary_{stem}.txt"

    np.savez_compressed(out_npz, **data)
    write_summary(out_txt, h5_path, timestep, data)

    print(f"Saved {out_npz}")
    print(f"Saved {out_txt}")


if __name__ == "__main__":
    main()
