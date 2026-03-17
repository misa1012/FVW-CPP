#!/usr/bin/env python3
"""
Project extracted filament segments to a 3D Cartesian grid as a regularized vorticity field.

Input:
    segments_timestep_*.npz produced by extract_filament_segments.py

Method:
    - Each filament segment contributes a vector strength:
          gamma * tangent * length
    - The segment center is deposited to the 8 surrounding grid nodes using
      trilinear (CIC-like) weights.
    - The nodal vector content is divided by a representative cell volume
      dx*dy*dz to obtain an equivalent vorticity field [1/s].

Outputs:
    - VTK image data (.vti) with point-data fields:
        omega
        omega_x
        omega_y
        omega_z
        omega_mag
    - summary text file
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter


TYPE_LABELS = {
    0: "bound",
    1: "shed",
    2: "trailing",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Project filament segments to a 3D vorticity grid")
    parser.add_argument("segments_npz", help="Path to segments_timestep_*.npz")
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory. Default: same directory as input npz",
    )
    parser.add_argument("--nx", type=int, default=160, help="Number of grid nodes in x")
    parser.add_argument("--ny", type=int, default=96, help="Number of grid nodes in y")
    parser.add_argument("--nz", type=int, default=96, help="Number of grid nodes in z")
    parser.add_argument("--pad-x", type=float, default=0.5, help="Padding added to x-bounds [m]")
    parser.add_argument("--pad-y", type=float, default=0.25, help="Padding added to y-bounds [m]")
    parser.add_argument("--pad-z", type=float, default=0.25, help="Padding added to z-bounds [m]")
    parser.add_argument("--diameter", type=float, default=None, help="Rotor diameter D [m] for normalized bounds")
    parser.add_argument("--hub-z", type=float, default=0.8, help="Hub center z [m] for normalized z-bounds")
    parser.add_argument("--xlim-d", nargs=2, type=float, default=None, help="Override x bounds in x/D")
    parser.add_argument("--ylim-d", nargs=2, type=float, default=None, help="Override y bounds in y/D")
    parser.add_argument("--zlim-d", nargs=2, type=float, default=None, help="Override z bounds in (z-hub)/D")
    parser.add_argument("--ppd", type=float, default=None, help="Points per diameter. Overrides nx/ny/nz when normalized bounds are used")
    parser.add_argument(
        "--types",
        default="bound,shed,trailing",
        help="Comma-separated subset of line types to include: bound,shed,trailing",
    )
    parser.add_argument(
        "--subdivide-factor",
        type=float,
        default=0.0,
        help=(
            "If > 0, subdivide each segment so subsegment length <= "
            "subdivide_factor * min(grid spacing). Recommended 1.0-2.0."
        ),
    )
    parser.add_argument(
        "--kernel",
        choices=["cic", "gaussian", "gaussian_fast"],
        default="cic",
        help="Deposition kernel type",
    )
    parser.add_argument(
        "--gaussian-sigma-cells",
        type=float,
        default=1.0,
        help="Gaussian kernel sigma in units of grid spacing (used when --kernel gaussian)",
    )
    return parser.parse_args()


def load_segments(npz_path: Path, allowed_types: set[int]) -> dict[str, np.ndarray]:
    data = np.load(npz_path)
    line_type = data["line_type"].astype(int)
    mask = np.isin(line_type, list(allowed_types))
    if not np.any(mask):
        raise RuntimeError("No segments left after type filtering")
    return {k: data[k][mask] for k in data.files}


def type_string_to_ids(spec: str) -> set[int]:
    lookup = {v: k for k, v in TYPE_LABELS.items()}
    ids = set()
    for item in spec.split(","):
        key = item.strip().lower()
        if not key:
            continue
        if key not in lookup:
            raise ValueError(f"Unknown line type '{item}'. Use bound,shed,trailing")
        ids.add(lookup[key])
    if not ids:
        raise ValueError("No valid line types selected")
    return ids


def build_grid(center_xyz: np.ndarray, nx: int, ny: int, nz: int, pad_x: float, pad_y: float, pad_z: float):
    xyz_min = center_xyz.min(axis=0)
    xyz_max = center_xyz.max(axis=0)
    origin = np.array([xyz_min[0] - pad_x, xyz_min[1] - pad_y, xyz_min[2] - pad_z], dtype=float)
    upper = np.array([xyz_max[0] + pad_x, xyz_max[1] + pad_y, xyz_max[2] + pad_z], dtype=float)
    dims = np.array([nx, ny, nz], dtype=int)
    spacing = (upper - origin) / (dims - 1)
    return origin, spacing, dims, upper


def build_grid_from_normalized_bounds(
    D: float,
    hub_z: float,
    xlim_d,
    ylim_d,
    zlim_d,
    nx: int,
    ny: int,
    nz: int,
    ppd: float | None,
):
    origin = np.array([xlim_d[0] * D, ylim_d[0] * D, hub_z + zlim_d[0] * D], dtype=float)
    upper = np.array([xlim_d[1] * D, ylim_d[1] * D, hub_z + zlim_d[1] * D], dtype=float)
    if ppd is not None:
        nx = int(round((xlim_d[1] - xlim_d[0]) * ppd)) + 1
        ny = int(round((ylim_d[1] - ylim_d[0]) * ppd)) + 1
        nz = int(round((zlim_d[1] - zlim_d[0]) * ppd)) + 1
    dims = np.array([nx, ny, nz], dtype=int)
    spacing = (upper - origin) / (dims - 1)
    return origin, spacing, dims, upper


def deposit_cic(center_xyz: np.ndarray, vec_strength: np.ndarray, origin: np.ndarray, spacing: np.ndarray, dims: np.ndarray):
    nx, ny, nz = map(int, dims)
    omega_content = np.zeros((nx, ny, nz, 3), dtype=np.float32)

    xi = (center_xyz - origin[None, :]) / spacing[None, :]
    base = np.floor(xi).astype(int)
    frac = xi - base

    valid = np.all(base >= 0, axis=1) & np.all(base < (dims - 1), axis=1)
    if not np.all(valid):
        center_xyz = center_xyz[valid]
        vec_strength = vec_strength[valid]
        base = base[valid]
        frac = frac[valid]

    for off_x in (0, 1):
        wx = (1.0 - frac[:, 0]) if off_x == 0 else frac[:, 0]
        ix = base[:, 0] + off_x
        for off_y in (0, 1):
            wy = (1.0 - frac[:, 1]) if off_y == 0 else frac[:, 1]
            iy = base[:, 1] + off_y
            for off_z in (0, 1):
                wz = (1.0 - frac[:, 2]) if off_z == 0 else frac[:, 2]
                iz = base[:, 2] + off_z
                w = wx * wy * wz
                np.add.at(omega_content, (ix, iy, iz, 0), vec_strength[:, 0] * w)
                np.add.at(omega_content, (ix, iy, iz, 1), vec_strength[:, 1] * w)
                np.add.at(omega_content, (ix, iy, iz, 2), vec_strength[:, 2] * w)

    return omega_content, int(valid.sum()), int((~valid).sum())


def deposit_gaussian(
    center_xyz: np.ndarray,
    vec_strength: np.ndarray,
    origin: np.ndarray,
    spacing: np.ndarray,
    dims: np.ndarray,
    sigma_cells: float,
):
    nx, ny, nz = map(int, dims)
    omega_content = np.zeros((nx, ny, nz, 3), dtype=np.float32)

    sigma_xyz = sigma_cells * spacing
    # compact support radius = 3 sigma
    radius_xyz = 3.0 * sigma_xyz

    used = 0
    dropped = 0

    for c, s in zip(center_xyz, vec_strength):
        lo = np.floor((c - radius_xyz - origin) / spacing).astype(int)
        hi = np.ceil((c + radius_xyz - origin) / spacing).astype(int)
        lo = np.maximum(lo, 0)
        hi = np.minimum(hi, dims - 1)
        if np.any(lo > hi):
            dropped += 1
            continue

        ix = np.arange(lo[0], hi[0] + 1)
        iy = np.arange(lo[1], hi[1] + 1)
        iz = np.arange(lo[2], hi[2] + 1)

        x = origin[0] + ix * spacing[0]
        y = origin[1] + iy * spacing[1]
        z = origin[2] + iz * spacing[2]

        wx = np.exp(-0.5 * ((x - c[0]) / sigma_xyz[0]) ** 2)
        wy = np.exp(-0.5 * ((y - c[1]) / sigma_xyz[1]) ** 2)
        wz = np.exp(-0.5 * ((z - c[2]) / sigma_xyz[2]) ** 2)

        w = (
            wx[:, None, None]
            * wy[None, :, None]
            * wz[None, None, :]
        ).astype(np.float32)
        wsum = float(w.sum())
        if wsum <= 0.0:
            dropped += 1
            continue
        w /= wsum

        omega_content[np.ix_(ix, iy, iz, [0])] += (w[..., None] * s[0]).astype(np.float32)
        omega_content[np.ix_(ix, iy, iz, [1])] += (w[..., None] * s[1]).astype(np.float32)
        omega_content[np.ix_(ix, iy, iz, [2])] += (w[..., None] * s[2]).astype(np.float32)
        used += 1

    return omega_content, used, dropped


def subdivide_segments(
    start_xyz: np.ndarray,
    end_xyz: np.ndarray,
    vec_strength: np.ndarray,
    spacing: np.ndarray,
    subdivide_factor: float,
):
    if subdivide_factor <= 0.0:
        return 0.5 * (start_xyz + end_xyz), vec_strength

    h = float(np.min(spacing))
    target = subdivide_factor * h
    if target <= 0.0:
        return 0.5 * (start_xyz + end_xyz), vec_strength

    seg_vec = end_xyz - start_xyz
    seg_len = np.linalg.norm(seg_vec, axis=1)
    nsub = np.maximum(1, np.ceil(seg_len / target).astype(int))

    centers = []
    strengths = []
    for p0, p1, strength, ns in zip(start_xyz, end_xyz, vec_strength, nsub):
        if ns == 1:
            centers.append(0.5 * (p0 + p1))
            strengths.append(strength)
            continue
        d = (p1 - p0) / ns
        sub_strength = strength / ns
        for k in range(ns):
            c = p0 + (k + 0.5) * d
            centers.append(c)
            strengths.append(sub_strength)
    return np.asarray(centers, dtype=np.float32), np.asarray(strengths, dtype=np.float32)


def write_vti(out_vti: Path, omega: np.ndarray, origin: np.ndarray, spacing: np.ndarray, dims: np.ndarray):
    grid = pv.ImageData()
    grid.origin = tuple(origin.tolist())
    grid.spacing = tuple(spacing.tolist())
    grid.dimensions = tuple(int(v) for v in dims.tolist())

    omega_flat = omega.reshape(-1, 3, order="F")
    grid.point_data["omega"] = omega_flat
    grid.point_data["omega_x"] = omega_flat[:, 0]
    grid.point_data["omega_y"] = omega_flat[:, 1]
    grid.point_data["omega_z"] = omega_flat[:, 2]
    grid.point_data["omega_mag"] = np.linalg.norm(omega_flat, axis=1)
    grid.save(out_vti)


def write_summary(
    out_txt: Path,
    npz_path: Path,
    allowed_types: set[int],
    n_segments: int,
    used_segments: int,
    dropped_segments: int,
    origin: np.ndarray,
    upper: np.ndarray,
    spacing: np.ndarray,
    dims: np.ndarray,
    vec_strength: np.ndarray,
    omega: np.ndarray,
) -> None:
    omega_flat = omega.reshape(-1, 3)
    omega_mag = np.linalg.norm(omega_flat, axis=1)
    total_in = vec_strength.sum(axis=0)
    dvol = float(np.prod(spacing))
    total_out = omega_flat.sum(axis=0) * dvol

    type_names = ", ".join(TYPE_LABELS[t] for t in sorted(allowed_types))
    lines = [
        "3D circulation-to-grid projection summary",
        f"segments_npz: {npz_path}",
        f"line_types: {type_names}",
        "",
        f"n_segments_input: {n_segments}",
        f"n_segments_used: {used_segments}",
        f"n_segments_dropped_outside_grid: {dropped_segments}",
        "",
        f"grid_dims: {dims[0]} x {dims[1]} x {dims[2]}",
        f"origin: ({origin[0]:.6f}, {origin[1]:.6f}, {origin[2]:.6f})",
        f"upper: ({upper[0]:.6f}, {upper[1]:.6f}, {upper[2]:.6f})",
        f"spacing: ({spacing[0]:.6f}, {spacing[1]:.6f}, {spacing[2]:.6f})",
        f"cell_volume: {dvol:.6e}",
        "",
        "Integrated vector content check:",
        f"  input  sum(gamma*t*ds): ({total_in[0]:.6e}, {total_in[1]:.6e}, {total_in[2]:.6e})",
        f"  output sum(omega*dV):   ({total_out[0]:.6e}, {total_out[1]:.6e}, {total_out[2]:.6e})",
        "",
        f"omega_mag min: {omega_mag.min():.6e}",
        f"omega_mag max: {omega_mag.max():.6e}",
        f"omega_mag mean: {omega_mag.mean():.6e}",
        "",
        f"omega_x range: [{omega_flat[:,0].min():.6e}, {omega_flat[:,0].max():.6e}]",
        f"omega_y range: [{omega_flat[:,1].min():.6e}, {omega_flat[:,1].max():.6e}]",
        f"omega_z range: [{omega_flat[:,2].min():.6e}, {omega_flat[:,2].max():.6e}]",
    ]
    out_txt.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    npz_path = Path(args.segments_npz).resolve()
    if not npz_path.exists():
        raise FileNotFoundError(f"Missing NPZ: {npz_path}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else npz_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    allowed_types = type_string_to_ids(args.types)
    seg = load_segments(npz_path, allowed_types)

    center_xyz = seg["center_xyz"]
    gamma = seg["gamma"]
    tangent = seg["tangent_xyz"]
    length = seg["length"]
    start_xyz = seg["start_xyz"]
    end_xyz = seg["end_xyz"]

    vec_strength = gamma[:, None] * tangent * length[:, None]

    if args.xlim_d is not None or args.ylim_d is not None or args.zlim_d is not None:
        if args.diameter is None:
            raise ValueError("--diameter is required when using normalized bounds")
        xlim_d = args.xlim_d if args.xlim_d is not None else [float(center_xyz[:, 0].min() / args.diameter), float(center_xyz[:, 0].max() / args.diameter)]
        ylim_d = args.ylim_d if args.ylim_d is not None else [float(center_xyz[:, 1].min() / args.diameter), float(center_xyz[:, 1].max() / args.diameter)]
        zrel = (center_xyz[:, 2] - args.hub_z) / args.diameter
        zlim_d = args.zlim_d if args.zlim_d is not None else [float(zrel.min()), float(zrel.max())]
        origin, spacing, dims, upper = build_grid_from_normalized_bounds(
            D=args.diameter,
            hub_z=args.hub_z,
            xlim_d=xlim_d,
            ylim_d=ylim_d,
            zlim_d=zlim_d,
            nx=args.nx,
            ny=args.ny,
            nz=args.nz,
            ppd=args.ppd,
        )
    else:
        origin, spacing, dims, upper = build_grid(
            center_xyz,
            nx=args.nx,
            ny=args.ny,
            nz=args.nz,
            pad_x=args.pad_x,
            pad_y=args.pad_y,
            pad_z=args.pad_z,
        )

    deposit_centers, deposit_strength = subdivide_segments(
        start_xyz,
        end_xyz,
        vec_strength,
        spacing,
        args.subdivide_factor,
    )

    if args.kernel == "cic":
        omega_content, used_segments, dropped_segments = deposit_cic(
            deposit_centers,
            deposit_strength,
            origin,
            spacing,
            dims,
        )
    elif args.kernel == "gaussian":
        omega_content, used_segments, dropped_segments = deposit_gaussian(
            deposit_centers,
            deposit_strength,
            origin,
            spacing,
            dims,
            sigma_cells=args.gaussian_sigma_cells,
        )
    else:
        omega_content, used_segments, dropped_segments = deposit_cic(
            deposit_centers,
            deposit_strength,
            origin,
            spacing,
            dims,
        )

    dvol = np.prod(spacing)
    omega = (omega_content / dvol).astype(np.float32, copy=False)
    if args.kernel == "gaussian_fast":
        sigma = float(args.gaussian_sigma_cells)
        for c in range(3):
            omega[..., c] = gaussian_filter(omega[..., c], sigma=sigma, mode="nearest")

    stem = npz_path.stem.replace("segments_", f"omega_grid_{args.kernel}_")
    out_vti = out_dir / f"{stem}.vti"
    out_txt = out_dir / f"{stem}_summary.txt"

    write_vti(out_vti, omega, origin, spacing, dims)
    write_summary(
        out_txt,
        npz_path,
        allowed_types,
        n_segments=len(length),
        used_segments=used_segments,
        dropped_segments=dropped_segments,
        origin=origin,
        upper=upper,
        spacing=spacing,
        dims=dims,
        vec_strength=vec_strength,
        omega=omega,
    )

    print(f"Saved {out_vti}")
    print(f"Saved {out_txt}")


if __name__ == "__main__":
    main()
