"""V8 stop-position sweep for an existing 5-surface V7 lens with a separate STOP.

Run from the project root. Requires optimize_multifield_v6.py in src/.
The input .zos MUST already contain:
  Surface 1: planar dummy STOP (air), with thickness to front lens
  Surface 2: front Even Asphere (N-BK7)
  Surface 3: rear lens surface
  Surface 4: image
No original input or earlier V7/V6 files are overwritten.
"""
import argparse
import csv
import math
from pathlib import Path

from optimize_multifield_v6 import connect, DEFAULT_ROOT, get_rms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='Saved .zos model with a separate dummy STOP at surface 1')
    parser.add_argument('--outdir', default='results/v8_stop_sweep')
    parser.add_argument('--zemax-root', default=DEFAULT_ROOT)
    parser.add_argument('--cfg1', type=Path, required=True)
    parser.add_argument('--cfg2', type=Path, required=True)
    parser.add_argument('--cfg3', type=Path, required=True)
    parser.add_argument('--start', type=float, default=0.0)
    parser.add_argument('--end', type=float, default=15.0)
    parser.add_argument('--step', type=float, default=1.0)
    parser.add_argument('--rms0-max', type=float, default=0.07)
    parser.add_argument('--rms2-max', type=float, default=0.238)
    args = parser.parse_args()
    src = Path(args.input).resolve()
    if not src.is_file():
        parser.error(f'Input file not found: {src}')
    cfg = {i: p.resolve() for i, p in enumerate((args.cfg1, args.cfg2, args.cfg3), 1)}
    for i, p in cfg.items():
        if not p.is_file():
            parser.error(f'CFG for field {i} not found: {p}')
    if args.start < 0 or args.step <= 0 or args.end < args.start:
        parser.error('Require 0 <= start <= end and step > 0')
    count = int(math.floor((args.end - args.start) / args.step + 1e-9)) + 1
    if count > 201:
        parser.error('At most 201 positions per run')
    out = Path(args.outdir).resolve()
    if out == src.parent:
        parser.error('Choose an output directory different from the input directory')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'reports').mkdir(exist_ok=True)
    app = None
    try:
        app, api = connect(args.zemax_root)
        system = app.PrimarySystem
        system.LoadFile(str(src), False)
        lde = system.LDE
        if lde.NumberOfSurfaces != 5:
            raise RuntimeError('Expected exactly 5 surfaces: OBJECT, STOP, Even Asphere, rear, IMAGE. Save the modified V7 lens first.')
        stop = lde.GetSurfaceAt(1)
        front = lde.GetSurfaceAt(2)
        rear = lde.GetSurfaceAt(3)
        # Do not silently sweep the thickness of an optical surface.
        if not math.isinf(float(stop.Radius)) or str(stop.Material).strip():
            raise RuntimeError('Surface 1 must be a plane dummy STOP with no material.')
        if str(front.Material).strip().upper() != 'N-BK7':
            raise RuntimeError('Expected N-BK7 front lens at surface 2. Check the LDE before running.')
        try:
            actual_stop = int(lde.StopSurface)
        except AttributeError:
            raise RuntimeError('Could not verify STOP surface via LDE.StopSurface; aborting for safety.')
        if actual_stop != 1:
            raise RuntimeError(f'STOP must be surface 1; current STOP is surface {actual_stop}.')
        original = (float(front.Radius), float(front.Thickness), float(rear.Radius), float(rear.Thickness))
        rows = []
        best = None
        print(f'Loaded: {src}\nVerified STOP=1; front lens=2; rear lens=3', flush=True)
        for i in range(count):
            d = args.start + i * args.step
            stop.Thickness = d
            rms = [float(get_rms(system, api, k, out / 'reports' / f'stop_{i:03d}_field{k}.txt', cfg)) for k in (1, 2, 3)]
            # Semi-diameters should be auto-computed in the source lens; if they are
            # fixed manually, these columns reflect those fixed values instead.
            sd_front = float(front.SemiDiameter)
            sd_rear = float(rear.SemiDiameter)
            feasible = rms[0] <= args.rms0_max and rms[1] <= args.rms2_max
            row = [round(d, 9), *rms, sd_front, sd_rear, 2 * max(sd_front, sd_rear), feasible]
            rows.append(row)
            print(f'STOP={d:5.2f} mm  RMS=[{rms[0]:.6f}, {rms[1]:.6f}, {rms[2]:.6f}]  '
                  f'diameter={2 * max(sd_front, sd_rear):.3f} mm  feasible={feasible}', flush=True)
            if feasible and (best is None or rms[2] < best[1]):
                best = (d, rms[2])
            if any(abs(v - w) > 1e-8 for v, w in zip(original, (float(front.Radius), float(front.Thickness), float(rear.Radius), float(rear.Thickness)))):
                raise RuntimeError('Unexpected lens geometry change during sweep; stopped.')
        csv_path = out / 'stop_sweep.csv'
        with csv_path.open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['STOP_distance_mm', 'RMS_0deg_waves', 'RMS_2deg_waves', 'RMS_5deg_waves',
                             'front_semi_diameter_mm', 'rear_semi_diameter_mm',
                             'required_lens_diameter_mm', 'feasible_RMS0_RMS2'])
            writer.writerows(rows)
        if best is not None:
            stop.Thickness = best[0]
            system.SaveAs(str(out / 'V8_best_stop.zos'))
            print(f'Best feasible STOP={best[0]:.3f} mm, RMS5={best[1]:.6f} waves', flush=True)
        else:
            print('No position satisfied both RMS constraints.', flush=True)
        print(f'Saved: {csv_path}', flush=True)
        print('Reminder: check ray clearance/vignetting and EFL in OpticStudio before accepting the design.', flush=True)
    finally:
        if app is not None:
            app.CloseApplication()


if __name__ == '__main__':
    main()
