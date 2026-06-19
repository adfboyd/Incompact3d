#!/usr/bin/env pvbatch
"""
Batch-render Xcompact3d ellipsoid snapshots.

Usage:
    pvbatch paraview_render.py data/snapshot-*.xdmf --outdir frames

Edit the USER SETTINGS block below — nothing else should need changing.
"""

# ================================================================
#  USER SETTINGS
# ================================================================

# --- Output ---
RESOLUTION  = (2160, 2160)   # (width, height) in pixels

# --- Slice renders ---
SLICE_NORMALS  = ["x", "y", "z"]  # which planes to render; subset e.g. ["x"]
SLICE_COORD    = None             # pin slice position (None = auto: uses far-body centroid)

# Colour ranges for slice fields (None = per-frame auto)
VORT_RANGE  = None   # fixed range e.g. (-5.0, 5.0), or None for auto
UX_RANGE    = None   # fixed range e.g. (-1.0, 2.0), or None for symmetric auto

# --- Body surface ---
BODY_COLOUR = (1.0, 1.0, 1.0)   # RGB; always fully opaque

# --- Standalone iso-surface renders (body + semi-transparent surface, no slice) ---
STANDALONE_OPACITY      = 0.6
VORT_STANDALONE_LEVELS  = [1.0, 1.5, 2.0]   # vorticity contour thresholds
CRITQ_STANDALONE_LEVELS = [-1.0, -2.0]             # Q-criterion thresholds; [] to disable
VORT_ISO_COLOUR         = (0.4, 0.75, 1.0)  # colour for vorticity iso-surfaces
CRITQ_ISO_COLOUR        = (0.2, 0.8,  0.3)  # colour for Q-criterion iso-surfaces

# --- Slice-overlaid vorticity iso-surface ---
VORT_ISO_THRESHOLD = None   # set to a float e.g. 3.0 to enable; None = disabled
VORT_ISO_OPACITY   = 0.5

# --- Rotating-camera orbit frame ---
# One frame per snapshot is saved alongside the standard views; the angle
# advances across all snapshots.  Set ROTATE_TOTAL = 0 to disable.
ROTATE_TOTAL = 150    # total degrees swept; 360 = full orbit, 90 = quarter turn
ROTATE_AXES  = ["x","y","z"]  # one rotating frame per entry; e.g. ["x", "y", "z"]

# --- Cinematic orbit ---
# When True: adds a 4th rotating frame per surface using an ease-in/out sweep
# that arcs the camera from CINEMATIC_ELEV_START down to CINEMATIC_ELEV_END
# degrees above horizontal (crane-shot feel), orbiting around CINEMATIC_AXIS.
# The 3 flat axis orbits (ROTATE_AXES) are always produced regardless.
# Set to False to produce only the 3 flat orbit frames per surface.
CINEMATIC            = True
CINEMATIC_AXIS       = "y"    # which axis the cinematic orbit starts on
CINEMATIC_ELEV_START = 30.0   # elevation at start of orbit (degrees above horizontal)
CINEMATIC_ELEV_END   =  8.0   # elevation at end of orbit

# --- Camera ---
VIEW_ANGLE = 30.0   # perspective FOV in degrees (applies to all renders)

# ================================================================
#  END OF USER SETTINGS
# ================================================================

import argparse
import math
import os
import re
import sys
import time

from paraview.simple import (
    XDMFReader, Calculator, Contour, Connectivity, Slice, Show, Hide, ColorBy,
    CreateRenderView, GetColorTransferFunction, SaveScreenshot,
    GetActiveCamera, RedistributeDataSet,
)
from paraview import servermanager
import numpy as np
from vtk.util import numpy_support as ns


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("xdmf", nargs="+",
                   help="snapshot-*.xdmf files (or a single file with embedded timesteps)")
    p.add_argument("--outdir", default="frames",
                   help="output directory for PNG frames (default: frames)")
    p.add_argument("--rotate-offset", type=int, default=0,
                   help="global index of the first snapshot in this shard (default: 0)")
    p.add_argument("--rotate-stride", type=int, default=1,
                   help="step between consecutive local snapshots in global sequence (default: 1)")
    p.add_argument("--total-frames", type=int, default=None,
                   help="total snapshots across all shards; if omitted uses len(xdmf)")
    return p.parse_args()


def snapshot_index(stem):
    """Parse the trailing integer from a snapshot stem, e.g. 'snapshot-123' -> 123.

    Xcompact3d writes snapshot numbers with Fortran '(I0)' (no zero-padding), so a
    shell glob delivers them in *lexicographic* order, not numeric ('snapshot-10'
    sorts before 'snapshot-2'). Driving the camera angle off the parsed number
    instead of the file's position in the shard makes the orbit immune to that.
    Returns None if no trailing integer is found.
    """
    m = re.search(r'(\d+)$', stem)
    return int(m.group(1)) if m else None


def ordered_partition(src):
    """Wrap a source in RedistributeDataSet for valid IceT ordered compositing.

    When rendered across MPI ranks, any translucent geometry forces IceT into
    ordered (back-to-front) compositing, which needs a spatially-sortable, non-empty
    partition on every rank. Iso-surfaces are sparse, so a domain-decomposed run
    leaves many ranks with empty geometry and degenerate bounds -> IceT aborts with
    'Invalid composite order'. Redistributing the *extracted surface* into a balanced
    kd-tree partition fixes the ordering. The 10^9-point volume field is never
    redistributed -- only the small surface coming out of the contour -- so the cost
    is an all-to-all of a few million cells, negligible beside reading the snapshot.
    """
    return RedistributeDataSet(Input=src)


def body_centroids(contour):
    """Return per-component (x,y,z) centroids from an ep1=0.5 contour.

    Uses ConnectivityFilter (already attached upstream) -> Fetch -> numpy.
    Returns an (M, 3) numpy array; empty (0, 3) array if no bodies present.
    """
    contour.UpdatePipeline()
    data = servermanager.Fetch(contour)
    n_pts = data.GetNumberOfPoints()
    if n_pts == 0:
        return np.zeros((0, 3))

    pts = ns.vtk_to_numpy(data.GetPoints().GetData())
    region_arr = data.GetPointData().GetArray("RegionId")
    if region_arr is None:
        return pts.mean(axis=0, keepdims=True)

    region = ns.vtk_to_numpy(region_arr)
    rids = np.unique(region)
    out = np.empty((len(rids), 3))
    for k, rid in enumerate(rids):
        out[k] = pts[region == rid].mean(axis=0)
    return out


def log(msg):
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()


def fmt_hms(seconds):
    """Format a duration in seconds as H:MM:SS."""
    seconds = int(round(max(seconds, 0)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    log(f"[start] {len(args.xdmf)} xdmf input(s), outdir={args.outdir}")

    reader = XDMFReader(FileNames=[args.xdmf[0]])
    log("[reader] XDMFReader created")
    available = list(reader.PointArrayStatus)
    want = ["ep1", "vort", "ux", "uy", "uz"]
    if CRITQ_STANDALONE_LEVELS:
        want.append("critq")
    needed = [n for n in want if n in available]
    reader.PointArrayStatus = needed
    for required in ("ep1", "vort", "ux"):
        if required not in needed:
            sys.stderr.write(f"ERROR: snapshots must contain '{required}' point array "
                             f"(found: {available})\n")
            sys.exit(1)

    # Body surface (ep1 = 0.5 iso-contour)
    body = Contour(Input=reader)
    body.ContourBy = ["POINTS", "ep1"]
    body.Isosurfaces = [0.5]
    body.ComputeNormals = 1

    conn = Connectivity(Input=body)
    conn.ColorRegions = 1

    # Mask fluid fields so values inside bodies don't skew auto-range or colouring
    calc_vort = Calculator(Input=reader)
    calc_vort.AttributeType = 'Point Data'
    calc_vort.ResultArrayName = 'vort_masked'
    calc_vort.Function = '(1-ep1)*vort'

    calc_ux = Calculator(Input=calc_vort)
    calc_ux.AttributeType = 'Point Data'
    calc_ux.ResultArrayName = 'ux_masked'
    calc_ux.Function = '(1-ep1)*ux'

    view_up_map = {"x": [0.0, 1.0, 0.0],
                   "y": [0.0, 0.0, 1.0],
                   "z": [0.0, 1.0, 0.0]}
    axis_idx_map = {"x": 0, "y": 1, "z": 2}

    slc = Slice(Input=calc_ux)
    slc.SliceType = "Plane"

    view = CreateRenderView()
    view.ViewSize = list(RESOLUTION)
    view.Background = [1.0, 1.0, 1.0]

    body_disp = Show(body, view)
    ColorBy(body_disp, None)
    body_disp.AmbientColor  = list(BODY_COLOUR)
    body_disp.DiffuseColor  = list(BODY_COLOUR)
    body_disp.Ambient       = 0.15
    body_disp.Diffuse       = 0.85
    body_disp.Specular      = 0.4
    body_disp.SpecularPower = 40
    body_disp.SpecularColor = [1.0, 1.0, 1.0]
    body_disp.Opacity       = 1.0

    def _make_standalone(vtk_field, tag_prefix, level, colour):
        lstr = str(int(level)) if level == int(level) else str(level)
        surf = Contour(Input=reader)
        surf.ContourBy = ["POINTS", vtk_field]
        surf.Isosurfaces = [level]
        surf.ComputeNormals = 1
        surf = ordered_partition(surf)   # translucent -> needs ordered compositing
        disp = Show(surf, view)
        ColorBy(disp, None)
        disp.AmbientColor  = list(colour)
        disp.DiffuseColor  = list(colour)
        disp.Ambient       = 0.15
        disp.Diffuse       = 0.85
        disp.Specular      = 0.3
        disp.SpecularPower = 30
        disp.Opacity       = STANDALONE_OPACITY
        Hide(surf, view)
        return (f"{tag_prefix}{lstr}", surf)

    standalone_surfs = []
    for level in VORT_STANDALONE_LEVELS:
        standalone_surfs.append(_make_standalone("vort", "viso", level, VORT_ISO_COLOUR))
        log(f"[setup] standalone vort iso-surface at level={level}")
    if CRITQ_STANDALONE_LEVELS:
        if "critq" not in needed:
            log("[setup] WARNING: critq not in snapshot data — critq renders will be skipped")
        else:
            for level in CRITQ_STANDALONE_LEVELS:
                standalone_surfs.append(_make_standalone("critq", "critq", level, CRITQ_ISO_COLOUR))
                log(f"[setup] standalone critq iso-surface at level={level}")

    if VORT_ISO_THRESHOLD is not None:
        log(f"[vort-iso] creating slice-overlaid iso-surface at threshold={VORT_ISO_THRESHOLD}")
        vort_surf = Contour(Input=reader)
        vort_surf.ContourBy = ["POINTS", "vort"]
        vort_surf.Isosurfaces = [VORT_ISO_THRESHOLD]
        vort_surf.ComputeNormals = 1
        vort_surf = ordered_partition(vort_surf)   # translucent -> ordered compositing
        vort_surf_disp = Show(vort_surf, view)
        ColorBy(vort_surf_disp, None)
        vort_surf_disp.AmbientColor  = list(VORT_ISO_COLOUR)
        vort_surf_disp.DiffuseColor  = list(VORT_ISO_COLOUR)
        vort_surf_disp.Ambient       = 0.15
        vort_surf_disp.Diffuse       = 0.85
        vort_surf_disp.Specular      = 0.3
        vort_surf_disp.SpecularPower = 30
        vort_surf_disp.Opacity       = VORT_ISO_OPACITY

    slc_disp = Show(slc, view)
    vortLUT = GetColorTransferFunction("vort_masked")
    vortLUT.ApplyPreset("Cool to Warm", True)
    uxLUT = GetColorTransferFunction("ux_masked")
    uxLUT.ApplyPreset("Cool to Warm (Diverging)", True)

    log("[reader] loading first snapshot to get domain bounds (this is the big read)...")
    reader.UpdatePipeline()
    log("[reader] first snapshot loaded")
    bounds = reader.GetDataInformation().GetBounds()
    centre = [0.5 * (bounds[2*k] + bounds[2*k+1]) for k in range(3)]
    extent = [bounds[2*k+1] - bounds[2*k] for k in range(3)]
    view.CameraParallelProjection = 0
    cam = GetActiveCamera()
    cam.SetViewAngle(VIEW_ANGLE)
    width_px, height_px = RESOLUTION
    aspect = width_px / float(height_px)
    half_fov = math.radians(VIEW_ANGLE / 2.0)

    # Build camera config for every axis needed by slices, rotation, or cinematic.
    axis_cfgs = {}
    for ax in dict.fromkeys(SLICE_NORMALS + ROTATE_AXES + [CINEMATIC_AXIS]):   # union, order-preserving
        ai = axis_idx_map[ax]
        nvec = [0.0, 0.0, 0.0]; nvec[ai] = 1.0
        vu = view_up_map[ax]
        vu_idx = vu.index(1.0)
        horiz_idx = next(k for k in (0, 1, 2) if k != ai and k != vu_idx)
        half_scene = max(0.5 * extent[vu_idx], 0.5 * extent[horiz_idx] / aspect)
        cam_dist = half_scene / math.tan(half_fov)
        axis_cfgs[ax] = dict(idx=ai, normal=nvec, view_up=vu, cam_offset=cam_dist)

    nfiles = len(args.xdmf)
    t0 = time.time()
    for i, xdmf_path in enumerate(args.xdmf):
        if i > 0:
            reader.FileNames = [xdmf_path]
            reader.UpdatePipeline()

        stem = os.path.splitext(os.path.basename(xdmf_path))[0]

        cents = body_centroids(conn)
        if cents.shape[0] == 0:
            log(f"[{i+1}/{nfiles}] {stem}  no body found, skipping")
            continue

        # --- Slice renders ---
        for ax in SLICE_NORMALS:
            cfg = axis_cfgs[ax]
            ai = cfg["idx"]
            # Camera sits on the +axis side, so the far body has the minimum
            # coordinate on this axis — slicing through it keeps both bodies visible.
            far_idx = int(np.argmin(cents[:, ai]))
            coord = SLICE_COORD if SLICE_COORD is not None else float(cents[far_idx, ai])

            slc.SliceType.Normal = cfg["normal"]
            cam.SetViewUp(*cfg["view_up"])

            for body_label, coord in [("", coord)]:
                origin = list(centre); origin[ai] = coord
                slc.SliceType.Origin = origin
                focal = list(centre); focal[ai] = coord
                pos = list(focal);    pos[ai] += cfg["cam_offset"]
                cam.SetFocalPoint(*focal)
                cam.SetPosition(*pos)
                slc.UpdatePipeline()

                if VORT_RANGE is None:
                    info = slc.GetPointDataInformation().GetArray("vort_masked")
                    if info is not None:
                        lo, hi = info.GetRange()
                        if hi > lo:
                            vortLUT.RescaleTransferFunction(lo, hi)
                else:
                    vortLUT.RescaleTransferFunction(*VORT_RANGE)

                if UX_RANGE is None:
                    info = slc.GetPointDataInformation().GetArray("ux_masked")
                    if info is not None:
                        lo, hi = info.GetRange()
                        m = max(abs(lo), abs(hi))
                        if m > 0:
                            uxLUT.RescaleTransferFunction(-m, m)
                else:
                    uxLUT.RescaleTransferFunction(*UX_RANGE)

                tag = f"_{body_label}" if body_label else ""
                for fld, lut, label in (("vort_masked", vortLUT, "vort"),
                                        ("ux_masked",   uxLUT,   "ux")):
                    ColorBy(slc_disp, ("POINTS", fld))
                    slc_disp.LookupTable = lut
                    slc_disp.SetScalarBarVisibility(view, False)
                    if VORT_ISO_THRESHOLD is not None:
                        Hide(vort_surf, view)
                        out = os.path.join(args.outdir, f"{stem}_n{ax}_{label}.png")
                        SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                        log(f"[{i+1}/{nfiles}] {stem}  n{ax}={coord:.4f}  field={label}  -> {out}")
                        Show(vort_surf, view)
                        out = os.path.join(args.outdir, f"{stem}_n{ax}_{label}_viso.png")
                        SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                        log(f"[{i+1}/{nfiles}] {stem}  n{ax}={coord:.4f}  field={label}+viso  -> {out}")
                    else:
                        out = os.path.join(args.outdir, f"{stem}_n{ax}_{label}.png")
                        SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                        log(f"[{i+1}/{nfiles}] {stem}  n{ax}={coord:.4f}  field={label}  -> {out}")

        if VORT_ISO_THRESHOLD is not None:
            Hide(vort_surf, view)

        # --- Standalone iso-surface renders (opaque body + semi-transparent surface) ---
        Hide(slc, view)
        if VORT_ISO_THRESHOLD is not None:
            Hide(vort_surf, view)

        for tag, surf in standalone_surfs:
            Show(surf, view)
            surf.UpdatePipeline()

            for ax in SLICE_NORMALS:
                cfg = axis_cfgs[ax]
                ai = cfg["idx"]
                cam.SetViewUp(*cfg["view_up"])
                focal = list(centre)
                pos = list(centre); pos[ai] += cfg["cam_offset"]
                cam.SetFocalPoint(*focal)
                cam.SetPosition(*pos)
                out = os.path.join(args.outdir, f"{stem}_{tag}_view{ax}.png")
                SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                log(f"[{i+1}/{nfiles}] {stem}  {tag}  view={ax}  -> {out}")

            Hide(surf, view)

        # --- Rotating-camera orbit (per surface: 3 flat axes + 1 cinematic = 4 per surface) ---
        total = args.total_frames if args.total_frames is not None else nfiles
        if ROTATE_TOTAL != 0 and total > 1:
            # Drive the orbit angle off the snapshot's own number, not its position
            # in this shard's file list. Snapshots use unpadded '(I0)' numbering, so a
            # glob hands them to us (and to the submit script's shards) in lexicographic
            # order; using the parsed number keeps the camera in true time order and
            # immune to that. Fall back to shard offset/stride if the name doesn't parse.
            snap_n = snapshot_index(stem)
            if snap_n is not None:
                global_i = snap_n
            else:
                global_i = args.rotate_offset + i * args.rotate_stride
            # Classic enumeration numbers snapshots 0..total-1, so total-1 is the span.
            # Clamp to guard against off-by-one / non-contiguous numbering.
            t   = min(max(global_i / float(total - 1), 0.0), 1.0)
            t_c = 0.5 - 0.5 * math.cos(math.pi * t)
            az_flat = t   * ROTATE_TOTAL
            az_cin  = t_c * ROTATE_TOTAL
            el_cin  = CINEMATIC_ELEV_START + (CINEMATIC_ELEV_END - CINEMATIC_ELEV_START) * t_c
            for tag, surf in standalone_surfs:
                Show(surf, view)
                surf.UpdatePipeline()
                # flat linear orbit — one frame per axis
                for rot_ax in ROTATE_AXES:
                    rot_cfg = axis_cfgs[rot_ax]
                    rot_ai  = rot_cfg["idx"]
                    cam.SetViewUp(*rot_cfg["view_up"])
                    focal = list(centre)
                    pos = list(centre); pos[rot_ai] += rot_cfg["cam_offset"]
                    cam.SetFocalPoint(*focal)
                    cam.SetPosition(*pos)
                    cam.Azimuth(az_flat)
                    out = os.path.join(args.outdir,
                                       f"{stem}_{tag}_rotating_n{rot_ax}.png")
                    SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                    log(f"[{i+1}/{nfiles}] {stem}  {tag}  rotating_n{rot_ax} "
                        f"az={az_flat:.1f}°  -> {out}")
                # cinematic orbit — eased azimuth + elevation arc
                if CINEMATIC:
                    cin_cfg = axis_cfgs[CINEMATIC_AXIS]
                    cin_ai  = cin_cfg["idx"]
                    cam.SetViewUp(*cin_cfg["view_up"])
                    focal = list(centre)
                    pos = list(centre); pos[cin_ai] += cin_cfg["cam_offset"]
                    cam.SetFocalPoint(*focal)
                    cam.SetPosition(*pos)
                    cam.Azimuth(az_cin)
                    cam.Elevation(el_cin)
                    cam.OrthogonalizeViewUp()
                    out = os.path.join(args.outdir,
                                       f"{stem}_{tag}_rotating_cinematic.png")
                    SaveScreenshot(out, view, ImageResolution=list(RESOLUTION))
                    log(f"[{i+1}/{nfiles}] {stem}  {tag}  rotating_cinematic "
                        f"az={az_cin:.1f}° el={el_cin:.1f}°  -> {out}")
                Hide(surf, view)

        # Restore slice for next snapshot (body was never hidden).
        Show(slc, view)
        if VORT_ISO_THRESHOLD is not None:
            Show(vort_surf, view)

        # --- Per-shard progress / ETA ---
        done = i + 1
        elapsed = time.time() - t0
        rate = elapsed / done
        log(f"[progress] {done}/{nfiles} snapshots  elapsed={fmt_hms(elapsed)}  "
            f"avg={rate:.1f}s/snap  ETA={fmt_hms(rate * (nfiles - done))}  "
            f"(est. total {fmt_hms(rate * nfiles)})")


if __name__ == "__main__":
    main()
