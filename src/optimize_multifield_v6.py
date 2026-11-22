"""Multi-field even-asphere optimization with ZOS-API.

Requires build_lens.py (connect, DEFAULT_ROOT), scipy, and a V5 .zos file.
Uses Zernike Standard analysis settings via typed .NET interface or a saved
analysis CFG file. Never silently substitutes an incorrect field.
"""
import argparse
import csv
import math
import re
from pathlib import Path
from scipy.optimize import minimize
from build_lens import connect, DEFAULT_ROOT

N = 1.518522
T = 5.0
APERTURE_RADIUS = 5.0
FIELDS = (0.0, 2.0, 5.0)
RMS_PATTERN = re.compile(r'RMS\s*\(to\s+centroid\)\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*waves', re.I)
FIELD_PATTERN = re.compile(r'\bField\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*\(deg\)', re.I)

def decode(path):
    data = path.read_bytes()
    if data.startswith((b'\xff\xfe', b'\xfe\xff')):
        return data.decode('utf-16')
    if b'\x00' in data[:200]:
        return data.decode('utf-16-le')
    return data.decode('utf-8-sig', errors='replace')

def solve_r2(r1, target=50.0):
    p1 = (N - 1) / r1
    den = 1 - (T / N) * p1
    if abs(den) < 1e-12:
        raise ValueError('Degenerate rear-surface power')
    p2 = (1 / target - p1) / den
    if abs(p2) < 1e-12:
        raise ValueError('Rear surface would be planar')
    return (1 - N) / p2

def sag(r, radius, a4=0.0, a6=0.0):
    q = 1 - (r / radius) ** 2
    if q <= 0:
        return math.nan
    return (r * r / radius) / (1 + math.sqrt(q)) + a4 * r**4 + a6 * r**6

def slope(r, radius, a4=0.0, a6=0.0):
    q = 1 - (r / radius) ** 2
    if q <= 0:
        return math.nan
    return r / (radius * math.sqrt(q)) + 4 * a4 * r**3 + 6 * a6 * r**5

def set_field(analysis, settings, field_num, cfg_files):
    """Select a field using an exposed typed interface or user-exported CFG.

    Zernike settings are sometimes exposed only as IAS_ in pythonnet.
    Never assume generic settings.Field exists.
    """
    candidates = [settings]
    for attr in ('IAS_ZernikeStandardCoefficients', 'IAS_ZernikeStandard'):
        obj = getattr(settings, attr, None)
        if obj is not None:
            candidates.append(obj)
    for obj in candidates:
        selector = getattr(obj, 'Field', None)
        if selector is not None:
            setter = getattr(selector, 'SetFieldNumber', None)
            if callable(setter):
                setter(field_num)
                return 'typed Field.SetFieldNumber'
        # Some Zernike interfaces expose a scalar Field property.
        try:
            value = getattr(obj, 'Field', None)
            if isinstance(value, (int, float)):
                setattr(obj, 'Field', field_num)
                if int(getattr(obj, 'Field')) == field_num:
                    return 'typed scalar Field'
        except (AttributeError, TypeError, ValueError):
            pass
    cfg = cfg_files.get(field_num)
    if cfg is not None:
        if not cfg.is_file():
            raise RuntimeError(f'Missing field {field_num} CFG: {cfg}')
        loader = getattr(settings, 'LoadFrom', None)
        if not callable(loader):
            raise RuntimeError('Analysis settings do not expose LoadFrom for CFG files')
        status = loader(str(cfg))
        if status is False:
            raise RuntimeError(f'Failed to load Zernike CFG: {cfg}')
        return f'CFG {cfg.name}'
    typ = settings.GetType().FullName if hasattr(settings, 'GetType') else type(settings).__name__
    raise RuntimeError(
        f'Cannot select Zernike Field {field_num} using settings type {typ}. '
        'Export three Zernike settings CFG files from OpticStudio (one for '
        'each Field 1, 2, 3) and supply --cfg1, --cfg2, --cfg3. '
        'The program will not reuse the default field.'
    )

def get_rms(system, zosapi, field_num, report, cfg_files):
    enum = zosapi.Analysis.AnalysisIDM
    kind = next((getattr(enum, name) for name in ('ZernikeStandardCoefficients', 'ZernikeStandard') if hasattr(enum, name)), None)
    if kind is None:
        raise RuntimeError('Zernike Standard analysis unavailable')
    analysis = system.Analyses.New_Analysis(kind)
    try:
        settings = analysis.GetSettings()
        method = set_field(analysis, settings, field_num, cfg_files)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.unlink(missing_ok=True)
        analysis.ApplyAndWaitForCompletion()
        analysis.GetResults().GetTextFile(str(report))
        if not report.exists() or report.stat().st_size == 0:
            raise RuntimeError(f'Empty Zernike report: {report}')
        txt = decode(report)
        match = RMS_PATTERN.search(txt)
        if not match:
            raise RuntimeError(f'RMS (to centroid) missing: {report}')
        fm = FIELD_PATTERN.search(txt)
        if not fm:
            raise RuntimeError(f'Cannot verify field angle in report: {report}')
        actual = float(fm.group(1))
        expected = FIELDS[field_num - 1]
        if abs(actual - expected) > 0.01:
            raise RuntimeError(f'Wrong field: expected {expected} deg, report says {actual} deg ({method})')
        value = float(match.group(1))
        if not math.isfinite(value) or value < 0:
            raise RuntimeError('Invalid wavefront RMS')
        return value
    finally:
        analysis.Close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--zemax-root', default=DEFAULT_ROOT)
    parser.add_argument('--input', default='results/v5/optimized_asphere.zos')
    parser.add_argument('--outdir', default='results/v6')
    parser.add_argument('--cfg1', type=Path)
    parser.add_argument('--cfg2', type=Path)
    parser.add_argument('--cfg3', type=Path)
    parser.add_argument('--maxiter', type=int, default=45)
    parser.add_argument('--min-edge', type=float, default=1.0)
    parser.add_argument('--max-slope', type=float, default=0.5)
    args = parser.parse_args()
    src = Path(args.input).resolve()
    if not src.is_file():
        parser.error(f'Input V5 model not found: {src}')
    cfg_files = {i: p.resolve() for i, p in enumerate((args.cfg1, args.cfg2, args.cfg3), 1) if p is not None}
    if cfg_files and len(cfg_files) != 3:
        parser.error('Supply all three of --cfg1, --cfg2, --cfg3, or none')
    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    app = None
    rows = []
    try:
        app, zosapi = connect(args.zemax_root)
        system = app.PrimarySystem
        system.LoadFile(str(src), False)
        fields = system.SystemData.Fields
        fields.SetFieldType(zosapi.SystemData.FieldType.Angle)
        while fields.NumberOfFields < 3:
            fields.AddField(0.0, 0.0, 1.0)
        while fields.NumberOfFields > 3:
            fields.DeleteField(fields.NumberOfFields)
        for i, angle in enumerate(FIELDS, 1):
            f = fields.GetField(i)
            f.X, f.Y, f.Weight = 0.0, angle, 1.0
        front = system.LDE.GetSurfaceAt(1)
        rear = system.LDE.GetSurfaceAt(2)
        col = zosapi.Editors.LDE.SurfaceColumn
        a4cell = front.GetSurfaceCell(col.Par2)
        a6cell = front.GetSurfaceCell(col.Par3)
        r1_initial = float(front.Radius)
        a4_initial = float(a4cell.DoubleValue)
        a6_initial = float(a6cell.DoubleValue)
        bfl_initial = float(rear.Thickness)
        x0 = [r1_initial, a4_initial / 1e-6, a6_initial / 1e-9, bfl_initial]
        bounds = [(25.0, 110.0), (-20.0, 20.0), (-10.0, 10.0), (35.0, 110.0)]
        def evaluate(x, tag, weights, record=True):
            r1, a4scaled, a6scaled, bfl = map(float, x)
            a4, a6 = a4scaled * 1e-6, a6scaled * 1e-9
            if any(not math.isfinite(v) or v < lo or v > hi for v, (lo, hi) in zip(x, bounds)):
                return 1e3
            try:
                r2 = solve_r2(r1)
                if not (math.isfinite(r2) and abs(r2) > APERTURE_RADIUS):
                    return 1e3
                radii = [APERTURE_RADIUS * k / 20 for k in range(21)]
                edges = [T - sag(r, r1, a4, a6) + sag(r, r2) for r in radii]
                slopes = [max(abs(slope(r, r1, a4, a6)), abs(slope(r, r2))) for r in radii]
                if not all(map(math.isfinite, edges + slopes)) or min(edges) < args.min_edge or max(slopes) > args.max_slope:
                    return 1e3
                front.Radius, rear.Radius = r1, r2
                a4cell.DoubleValue, a6cell.DoubleValue = a4, a6
                rear.Thickness = bfl
                idx = len(rows)
                rms = [get_rms(system, zosapi, k, out / 'reports' / f'{tag}_{idx:05d}_field{k}.txt', cfg_files) for k in (1, 2, 3)]
                score = math.sqrt(sum(w * v*v for w, v in zip(weights, rms)) / sum(weights))
                if record:
                    rows.append([tag, idx, r1, r2, a4, a6, bfl, min(edges), max(slopes), *rms, score])
                    if idx % 10 == 0:
                        print(f'{tag}: evaluation {idx}, RMS={rms}, objective={score:.6f}', flush=True)
                return score
            except Exception as exc:
                raise RuntimeError(f'ZOS-API evaluation failed in {tag}: {exc}') from exc
        base = evaluate(x0, 'V5_reference', (1, 1, 1))
        if base >= 1e3:
            raise RuntimeError('V5 baseline violates optimization bounds/geometry; adjust bounds before proceeding')
        print('V5 baseline weighted RMS:', base)
        system.SaveAs(str(out / 'v5_three_fields.zos'))
        for tag, weights in [('V6_A_balanced', (1, 1, 1)), ('V6_B_edge_priority', (1, 1, 4))]:
            best = {'score': math.inf, 'x': None}
            def objective(x):
                score = evaluate(x, tag, weights)
                if score < best['score']:
                    best.update(score=score, x=list(x))
                return score
            result = minimize(objective, x0, method='Powell', bounds=bounds,
                              options={'maxiter': args.maxiter, 'maxfev': 350, 'xtol': 1e-3, 'ftol': 1e-4})
            if best['x'] is None or best['score'] >= 1e3:
                raise RuntimeError(f'No valid candidate for {tag}')
            evaluate(best['x'], tag + '_FINAL', weights)
            system.SaveAs(str(out / (tag + '.zos')))
            print(f'{tag}: best weighted RMS={best["score"]:.6f}, success={result.success}, {result.message}')
        with (out / 'comparison.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['stage', 'evaluation', 'R1_mm', 'R2_mm', 'A4_mm-3', 'A6_mm-5', 'BFL_mm', 'minimum_thickness_mm', 'maximum_slope', 'RMS_0deg_waves', 'RMS_2deg_waves', 'RMS_5deg_waves', 'weighted_RMS_waves'])
            writer.writerows(rows)
        print('Saved V6 models and comparison.csv to', out)
        print('VERIFY: Zemax Prescription EFL, field reports, FFT MTF, wavefront, sag, and tolerances.')
    finally:
        if app is not None:
            app.CloseApplication()

if __name__ == '__main__':
    main()
