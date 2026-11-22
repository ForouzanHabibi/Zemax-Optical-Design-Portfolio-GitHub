"""V7 short constrained multi-field feasibility trial.

Place beside optimize_multifield_v6.py and build_lens.py in src/.
Uses V6's verified Zernike CFG/field-checking implementation.
No V5/V6 files are overwritten. No false claim of convergence.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from optimize_multifield_v6 import (
    connect, DEFAULT_ROOT, N, T, APERTURE_RADIUS, FIELDS,
    solve_r2, sag, slope, get_rms,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='results/v6/v5_three_fields.zos')
    p.add_argument('--outdir', default='results/v7')
    p.add_argument('--zemax-root', default=DEFAULT_ROOT)
    p.add_argument('--cfg1', type=Path)
    p.add_argument('--cfg2', type=Path)
    p.add_argument('--cfg3', type=Path)
    p.add_argument('--maxfev', type=int, default=100)
    p.add_argument('--rms0-max', type=float, default=0.07)
    p.add_argument('--rms2-max', type=float, default=0.238)
    p.add_argument('--min-edge', type=float, default=1.0)
    p.add_argument('--max-slope', type=float, default=0.5)
    args = p.parse_args()
    src = Path(args.input).resolve()
    if not src.is_file():
        p.error(f'Missing input: {src}')
    cfg = {i: q.resolve() for i, q in enumerate((args.cfg1, args.cfg2, args.cfg3), 1) if q is not None}
    if cfg and len(cfg) != 3:
        p.error('Supply all three CFG paths, or none')
    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    app = None
    try:
        app, api = connect(args.zemax_root)
        sys = app.PrimarySystem
        sys.LoadFile(str(src), False)
        fs = sys.SystemData.Fields
        fs.SetFieldType(api.SystemData.FieldType.Angle)
        while fs.NumberOfFields < 3:
            fs.AddField(0.0, 0.0, 1.0)
        while fs.NumberOfFields > 3:
            fs.DeleteField(fs.NumberOfFields)
        for k, angle in enumerate(FIELDS, 1):
            f = fs.GetField(k)
            f.X, f.Y, f.Weight = 0.0, angle, 1.0
        front, rear = sys.LDE.GetSurfaceAt(1), sys.LDE.GetSurfaceAt(2)
        col = api.Editors.LDE.SurfaceColumn
        c4 = front.GetSurfaceCell(col.Par2)
        c6 = front.GetSurfaceCell(col.Par3)
        x0 = np.array([float(front.Radius), float(c4.DoubleValue) / 1e-6,
                       float(c6.DoubleValue) / 1e-9, float(rear.Thickness)])
        bounds = [(25., 110.), (-20., 20.), (-10., 10.), (35., 110.)]
        rows = []
        best = {'score': math.inf, 'x': None, 'rms': None}
        best_feasible = {'rms5': math.inf, 'x': None, 'rms': None}
        cache = {}

        def apply_geometry(x):
            r1, a4s, a6s, bfl = [float(v) for v in x]
            if any(not math.isfinite(v) or v < lo or v > hi for v, (lo, hi) in zip(x, bounds)):
                return None
            a4, a6 = a4s * 1e-6, a6s * 1e-9
            try:
                r2 = solve_r2(r1)
                if not math.isfinite(r2) or abs(r2) <= APERTURE_RADIUS:
                    return None
                radii = [APERTURE_RADIUS * k / 20 for k in range(21)]
                edges = [T - sag(r, r1, a4, a6) + sag(r, r2) for r in radii]
                slopes = [max(abs(slope(r, r1, a4, a6)), abs(slope(r, r2))) for r in radii]
                if (not all(map(math.isfinite, edges + slopes)) or
                        min(edges) < args.min_edge or max(slopes) > args.max_slope):
                    return None
            except (ValueError, ZeroDivisionError, OverflowError):
                return None
            front.Radius, rear.Radius = r1, r2
            c4.DoubleValue, c6.DoubleValue = a4, a6
            rear.Thickness = bfl
            return (r1, r2, a4, a6, bfl, min(edges), max(slopes))

        def evaluate(x):
            key = tuple(round(float(v), 12) for v in x)
            if key in cache:
                return cache[key]
            geo = apply_geometry(x)
            if geo is None:
                return 1e4
            i = len(rows)
            rms = [get_rms(sys, api, k, out / 'reports' / f'v7_{i:04d}_field{k}.txt', cfg)
                   for k in (1, 2, 3)]
            # Quadratic exterior penalty: trial optimization, NOT proof of feasibility.
            violation0 = max(0., rms[0] - args.rms0_max)
            violation2 = max(0., rms[1] - args.rms2_max)
            score = rms[2] + 100. * (violation0**2 + violation2**2)
            feasible = violation0 <= 1e-8 and violation2 <= 1e-8
            rows.append([i, *geo, *rms, violation0, violation2, score, feasible])
            if score < best['score']:
                best.update(score=score, x=list(x), rms=rms)
            if feasible and rms[2] < best_feasible['rms5']:
                best_feasible.update(rms5=rms[2], x=list(x), rms=rms)
            if i % 10 == 0:
                print(f'V7 eval {i}: RMS={rms}, penalty_objective={score:.6f}, feasible={feasible}', flush=True)
            cache[key] = score
            return score

        base_score = evaluate(x0)
        print('V5 starting RMS:', rows[0][8:11], 'score:', base_score, flush=True)
        result = minimize(evaluate, x0, method='Powell', bounds=bounds,
                          options={'maxfev': args.maxfev, 'maxiter': 15,
                                   'xtol': 1e-3, 'ftol': 1e-4})
        header = ['evaluation', 'R1_mm', 'R2_mm', 'A4_mm-3', 'A6_mm-5',
                  'BFL_mm', 'minimum_thickness_mm', 'maximum_slope',
                  'RMS_0deg_waves', 'RMS_2deg_waves', 'RMS_5deg_waves',
                  'violation_0', 'violation_2', 'penalty_objective', 'feasible']
        with (out / 'v7_trial.csv').open('w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        if best['x'] is not None:
            apply_geometry(best['x'])
            sys.SaveAs(str(out / 'V7_best_penalized.zos'))
        if best_feasible['x'] is not None:
            apply_geometry(best_feasible['x'])
            sys.SaveAs(str(out / 'V7_best_feasible.zos'))
            print('FEASIBLE candidate found, RMS:', best_feasible['rms'], flush=True)
        else:
            print('NO FEASIBLE candidate found within this short trial; this does NOT prove infeasibility.', flush=True)
        print('Optimizer success:', result.success, 'reason:', result.message, flush=True)
        print('Best penalized RMS:', best['rms'], 'objective:', best['score'], flush=True)
        print('Saved results to:', out, flush=True)
    finally:
        if app is not None:
            app.CloseApplication()


if __name__ == '__main__':
    main()
