"""V9 joint STOP + even-asphere optimization for the 5-surface V8 lens.

Run from project root; put this script in src/ alongside optimize_multifield_v6.py
and build_lens.py. Original ZOS and earlier results are never overwritten.

This is a constrained numerical trial, not a guarantee of a global optimum.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from optimize_multifield_v6 import (connect, DEFAULT_ROOT, N, T, FIELDS,
                                    solve_r2, sag, slope, get_rms)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', default='results/v8_stop_sweep_40_80/V8_best_stop.zos')
    p.add_argument('--outdir', default='results/v9_joint')
    p.add_argument('--zemax-root', default=DEFAULT_ROOT)
    for i in (1, 2, 3):
        p.add_argument(f'--cfg{i}', type=Path, required=True)
    p.add_argument('--maxfev', type=int, default=250)
    p.add_argument('--stop-min', type=float, default=20.)
    p.add_argument('--stop-max', type=float, default=80.)
    p.add_argument('--diameter-max', type=float, default=22.)
    p.add_argument('--rms0-max', type=float, default=0.07)
    p.add_argument('--rms2-max', type=float, default=0.238)
    p.add_argument('--min-edge', type=float, default=1.)
    p.add_argument('--max-slope', type=float, default=0.5)
    p.add_argument('--penalty', type=float, default=1000.)
    args = p.parse_args()
    src = Path(args.input).resolve()
    if not src.is_file():
        p.error(f'Input ZOS missing: {src}')
    cfg = {i: getattr(args, f'cfg{i}').resolve() for i in (1, 2, 3)}
    for k, path in cfg.items():
        if not path.is_file():
            p.error(f'Field {k} CFG missing: {path}')
    if not (0 <= args.stop_min < args.stop_max and args.diameter_max >= 10 and
            args.maxfev >= 10 and args.penalty > 0):
        p.error('Invalid optimization bounds/limits')
    out = Path(args.outdir).resolve()
    if out == src.parent or src == out / src.name:
        p.error('Output directory must differ from input directory')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'reports').mkdir(exist_ok=True)
    app = None
    try:
        app, api = connect(args.zemax_root)
        system = app.PrimarySystem
        system.LoadFile(str(src), False)
        lde = system.LDE
        if lde.NumberOfSurfaces != 5 or int(lde.StopSurface) != 1:
            raise RuntimeError('Expected OBJECT, STOP, Even Asphere, rear, IMAGE; STOP must be surface 1')
        stop, front, rear = (lde.GetSurfaceAt(i) for i in (1, 2, 3))
        if not math.isinf(float(stop.Radius)) or str(stop.Material).strip():
            raise RuntimeError('Surface 1 must be a planar, air STOP')
        if str(front.Material).strip().upper() != 'N-BK7':
            raise RuntimeError('Surface 2 must be the N-BK7 even asphere')
        if abs(float(front.Thickness) - T) > 1e-5:
            raise RuntimeError(f'Expected center thickness {T} mm; found {front.Thickness}')
        fs = system.SystemData.Fields
        if fs.NumberOfFields != 3:
            raise RuntimeError('Expected exactly 3 fields (0, 2, 5 degrees)')
        for k, expected in enumerate(FIELDS, 1):
            f = fs.GetField(k)
            if abs(float(f.Y) - expected) > 1e-5 or abs(float(f.X)) > 1e-5:
                raise RuntimeError(f'Field {k} does not match {expected} deg')
        col = api.Editors.LDE.SurfaceColumn
        c4 = front.GetSurfaceCell(col.Par2)
        c6 = front.GetSurfaceCell(col.Par3)
        x0 = np.array([float(front.Radius), float(c4.DoubleValue) / 1e-6,
                       float(c6.DoubleValue) / 1e-9, float(rear.Thickness),
                       float(stop.Thickness)])
        bounds = [(25., 110.), (-20., 20.), (-10., 10.), (35., 110.),
                  (args.stop_min, args.stop_max)]
        if any(not (lo <= v <= hi) for v, (lo, hi) in zip(x0, bounds)):
            raise RuntimeError(f'Starting geometry outside bounds: {x0.tolist()}')
        rows, cache = [], {}
        best_score, best_feasible = None, None

        def apply_geometry(x):
            try:
                r1, a4s, a6s, bfl, stopdist = [float(v) for v in x]
                if any(not math.isfinite(v) or v < lo or v > hi
                       for v, (lo, hi) in zip(x, bounds)):
                    return None
                r2 = solve_r2(r1)
                a4, a6 = a4s * 1e-6, a6s * 1e-9
                radii = np.linspace(0., args.diameter_max / 2., 41)
                edges = [T - sag(r, r1, a4, a6) + sag(r, r2) for r in radii]
                slopes = [max(abs(slope(r, r1, a4, a6)), abs(slope(r, r2)))
                          for r in radii]
                if (not all(map(math.isfinite, edges + slopes)) or
                    min(edges) < args.min_edge or max(slopes) > args.max_slope):
                    return None
                stop.Thickness = stopdist
                front.Radius = r1
                rear.Radius = r2
                c4.DoubleValue = a4
                c6.DoubleValue = a6
                rear.Thickness = bfl
                return [r1, r2, a4, a6, bfl, stopdist, min(edges), max(slopes)]
            except (ValueError, ZeroDivisionError, OverflowError):
                return None

        def evaluate(x):
            nonlocal best_score, best_feasible
            key = tuple(round(float(v), 11) for v in x)
            if key in cache:
                return cache[key]
            geo = apply_geometry(x)
            if geo is None:
                return 1e6
            idx = len(rows)
            rms = [float(get_rms(system, api, k,
                                out / 'reports' / f'v9_{idx:04d}_field{k}.txt', cfg))
                   for k in (1, 2, 3)]
            # Evaluate the lens size demanded by actual traced rays, not merely
            # the user-entered aperture or the nominal 10-mm entrance pupil.
            sd_front = float(front.SemiDiameter)
            sd_rear = float(rear.SemiDiameter)
            diameter = 2. * max(sd_front, sd_rear)
            if not math.isfinite(diameter):
                raise RuntimeError('Invalid lens diameter from Zemax')
            violations = [max(0., rms[0] - args.rms0_max),
                          max(0., rms[1] - args.rms2_max),
                          max(0., (diameter - args.diameter_max) / args.diameter_max)]
            score = rms[2] + args.penalty * sum(v*v for v in violations)
            feasible = all(v <= 1e-8 for v in violations)
            rows.append([idx, *geo, sd_front, sd_rear, diameter, *rms,
                         *violations, score, feasible])
            if best_score is None or score < best_score['score']:
                best_score = dict(x=list(x), rms=rms, diameter=diameter, score=score)
            if feasible and (best_feasible is None or rms[2] < best_feasible['rms'][2]):
                best_feasible = dict(x=list(x), rms=rms, diameter=diameter, score=score)
            if idx % 10 == 0:
                print(f'V9 eval {idx}: STOP={geo[5]:.2f}, RMS={rms}, '
                      f'diameter={diameter:.3f}, feasible={feasible}', flush=True)
            cache[key] = score
            return score

        print('Starting vector [R1, A4/1e-6, A6/1e-9, BFL, STOP]:', x0.tolist(), flush=True)
        evaluate(x0)
        result = minimize(evaluate, x0, method='Powell', bounds=bounds,
                          options={'maxfev': args.maxfev, 'maxiter': 25,
                                   'xtol': 1e-3, 'ftol': 1e-4})
        header = ['evaluation', 'R1_mm', 'R2_mm', 'A4_mm^-3', 'A6_mm^-5',
                  'BFL_mm', 'STOP_mm', 'minimum_edge_mm', 'max_surface_slope',
                  'front_semi_diameter_mm', 'rear_semi_diameter_mm',
                  'required_diameter_mm', 'RMS_0deg_waves', 'RMS_2deg_waves',
                  'RMS_5deg_waves', 'violation_RMS0', 'violation_RMS2',
                  'violation_diameter_relative', 'penalty_objective', 'feasible']
        with (out / 'v9_trials.csv').open('w', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(header)
            writer.writerows(rows)
        if best_score:
            apply_geometry(best_score['x'])
            system.SaveAs(str(out / 'V9_best_penalized.zos'))
        if best_feasible:
            apply_geometry(best_feasible['x'])
            system.SaveAs(str(out / 'V9_best_feasible.zos'))
            print('Best FEASIBLE:', best_feasible, flush=True)
        else:
            print('No candidate met all RMS and diameter limits.', flush=True)
        print('Best penalized:', best_score, flush=True)
        print('Optimizer:', result.success, result.message, flush=True)
        print('Output:', out, flush=True)
        print('IMPORTANT: verify EFL, ray clearance, vignetting, and MTF in OpticStudio.', flush=True)
    finally:
        if app is not None:
            app.CloseApplication()


if __name__ == '__main__':
    main()
