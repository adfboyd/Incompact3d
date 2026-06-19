#!/usr/bin/env python3
"""
Compute masked fluid kinetic energy from Incompact3d snapshot files.

Reads ux, uy, uz, ep1 binary snapshots from a run directory and writes
a CSV with (snapshot_num, itime, t, ke) for each snapshot.  ep1 marks
the fictitious solid-interior flow, so the fluid mask is (1 - ep1).

KE = 0.5 * sum((1 - ep1) * (ux^2 + uy^2 + uz^2)) * dx * dy * dz

File format: raw little-endian float64, C-order (x varies fastest in
Fortran = slowest in numpy), shape (nz, ny, nx) when reshaped — but
for a plain sum the ordering is irrelevant, so we work on the flat
1-D array in chunks.

Usage:
    python3 compute_ke.py /path/to/run [--input input.i3d] [--chunk 10000000] [--out ke.csv]

The run directory must contain:
    data/ux-N.bin  data/uy-N.bin  data/uz-N.bin  data/ep1-N.bin

Parameters are read from input.i3d in the run directory (or --input).
"""

import argparse
import csv
import glob
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np


# ---------------------------------------------------------------------------
# Input-file parser
# ---------------------------------------------------------------------------

def _val(text, key):
    """Return the value of 'key = <value>' from an i3d file as a string."""
    m = re.search(
        r'^\s*' + re.escape(key) + r'\s*=\s*([^\s,!&/]+)',
        text, re.MULTILINE | re.IGNORECASE,
    )
    if m is None:
        raise ValueError(f"Could not find '{key}' in input file")
    return m.group(1)


def parse_input(path):
    with open(path) as f:
        text = f.read()

    nx      = int(_val(text, 'nx'))
    ny      = int(_val(text, 'ny'))
    nz      = int(_val(text, 'nz'))
    xlx     = float(_val(text, 'xlx'))
    yly     = float(_val(text, 'yly'))
    zlz     = float(_val(text, 'zlz'))
    dt      = float(_val(text, 'dt'))
    ioutput = int(_val(text, 'ioutput'))

    # Periodic grids: dx = L / n  (no endpoint duplication)
    dx = xlx / nx
    dy = yly / ny
    dz = zlz / nz

    return dict(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz, dt=dt, ioutput=ioutput)


# ---------------------------------------------------------------------------
# Chunked KE computation
# ---------------------------------------------------------------------------

def ke_from_snapshot(data_dir, num, n_total, chunk, dtype='<f8'):
    """
    Read ux/uy/uz/ep1 for snapshot *num* in chunks and return (ke_sum, elapsed).
    Uses numpy.fromfile (C-level, no GIL contention on allocation) so this
    function is safe to call from multiple threads concurrently.
    """
    paths = [
        os.path.join(data_dir, f'ux-{num}.bin'),
        os.path.join(data_dir, f'uy-{num}.bin'),
        os.path.join(data_dir, f'uz-{num}.bin'),
        os.path.join(data_dir, f'ep1-{num}.bin'),
    ]
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing: {p}")

    ke_sum = np.float64(0.0)
    t0 = time.perf_counter()
    handles = [open(p, 'rb') for p in paths]
    try:
        while True:
            ux = np.fromfile(handles[0], dtype=dtype, count=chunk)
            if len(ux) == 0:
                break
            uy  = np.fromfile(handles[1], dtype=dtype, count=chunk)
            uz  = np.fromfile(handles[2], dtype=dtype, count=chunk)
            ep1 = np.fromfile(handles[3], dtype=dtype, count=chunk)
            ke_sum += np.dot(1.0 - ep1, ux * ux + uy * uy + uz * uz)
    finally:
        for h in handles:
            h.close()

    return ke_sum, time.perf_counter() - t0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', help='Path to run directory containing data/')
    ap.add_argument('--input', default=None,
                    help='Path to input.i3d (default: run_dir/input.i3d)')
    ap.add_argument('--chunk', type=int, default=50_000_000,
                    help='Elements per field per chunk (default 50M = ~400 MB/field)')
    ap.add_argument('--snap-workers', type=int, default=4,
                    help='Snapshots to process in parallel (default 4)')
    ap.add_argument('--out', default='ke.csv',
                    help='Output CSV file (default: ke.csv in run_dir)')
    args = ap.parse_args()

    run_dir  = os.path.abspath(args.run_dir)
    data_dir = os.path.join(run_dir, 'data')

    if not os.path.isdir(data_dir):
        sys.exit(f"Error: data/ directory not found in {run_dir}")

    input_path = args.input or os.path.join(run_dir, 'input.i3d')
    if not os.path.exists(input_path):
        sys.exit(f"Error: input file not found: {input_path}")

    params = parse_input(input_path)
    nx, ny, nz   = params['nx'], params['ny'], params['nz']
    dx, dy, dz   = params['dx'], params['dy'], params['dz']
    dt, ioutput  = params['dt'], params['ioutput']
    n_total      = nx * ny * nz
    cell_vol     = 0.5 * dx * dy * dz

    print(f"Grid: {nx} x {ny} x {nz}  ({n_total/1e6:.1f}M points, "
          f"{4 * n_total * 8 / 1e9:.1f} GB per snapshot)")
    print(f"dx={dx:.6g}  dy={dy:.6g}  dz={dz:.6g}  dt={dt:.6g}  ioutput={ioutput}")
    print(f"Chunk size: {args.chunk/1e6:.0f}M elements = "
          f"{4 * args.chunk * 8 / 1e6:.0f} MB in memory at once")

    # Find all snapshot numbers by globbing ep1 (always written, even for inviscid)
    ep1_files = sorted(
        glob.glob(os.path.join(data_dir, 'ep1-*.bin')),
        key=lambda p: int(re.search(r'ep1-(\d+)\.bin', p).group(1))
    )
    if not ep1_files:
        sys.exit("Error: no ep1-*.bin files found in data/")

    nums = [int(re.search(r'ep1-(\d+)\.bin', p).group(1)) for p in ep1_files]
    print(f"\nFound {len(nums)} snapshots: {nums[0]} ... {nums[-1]}\n")

    out_path = args.out if os.path.isabs(args.out) else os.path.join(run_dir, args.out)

    def fmt_duration(seconds):
        seconds = int(seconds)
        h, rem = divmod(seconds, 3600)
        m, s   = divmod(rem, 60)
        if h:
            return f"{h}h{m:02d}m{s:02d}s"
        if m:
            return f"{m}m{s:02d}s"
        return f"{s}s"

    n_snaps    = len(nums)
    snap_bytes = 4 * n_total * 8  # 4 fields × float64
    wall_start = time.perf_counter()

    # Progress state shared across threads
    progress_lock  = threading.Lock()
    completed      = {'n': 0}

    def print_progress(num, ke, elapsed):
        gbps = snap_bytes / elapsed / 1e9
        with progress_lock:
            completed['n'] += 1
            done         = completed['n']
            wall_elapsed = time.perf_counter() - wall_start
            avg          = wall_elapsed / done
            eta          = avg * (n_snaps - done)
            itime_       = num * ioutput
            t_           = itime_ * dt
            print(
                f"  [{done:4d}/{n_snaps}]  snap={num:6d}  t={t_:.4f}  KE={ke:.6e}"
                f"  {elapsed:5.1f}s  {gbps:.1f}GB/s"
                f"  wall={fmt_duration(wall_elapsed)}"
                f"  ETA={fmt_duration(eta)}",
                flush=True,
            )

    with ThreadPoolExecutor(max_workers=args.snap_workers) as executor, \
         open(out_path, 'w', newline='') as csvfile:

        writer = csv.writer(csvfile)
        writer.writerow(['snapshot_num', 'itime', 't', 'ke'])

        # Submit all snapshots; executor caps concurrent work at snap_workers
        future_to_num = {
            executor.submit(ke_from_snapshot, data_dir, num, n_total, args.chunk): num
            for num in nums
        }

        # Collect results as they finish (for progress), buffer for ordered CSV
        pending = {}
        next_idx = 0

        for fut in as_completed(future_to_num):
            num = future_to_num[fut]
            ke_sum, elapsed = fut.result()
            ke = ke_sum * cell_vol
            pending[num] = (ke, elapsed)
            print_progress(num, ke, elapsed)

            # Flush any consecutive completed rows to CSV in order
            while next_idx < n_snaps and nums[next_idx] in pending:
                n      = nums[next_idx]
                ke_row, _ = pending.pop(n)
                itime  = n * ioutput
                t      = itime * dt
                writer.writerow([n, itime, f'{t:.8g}', f'{ke_row:.12g}'])
                csvfile.flush()
                next_idx += 1

    total = time.perf_counter() - wall_start
    print(f"\nDone in {fmt_duration(total)} ({total/n_snaps:.1f}s/snap avg).  Written {out_path}", flush=True)


if __name__ == '__main__':
    main()
