"""V4: fixed-EFL spherical singlet shape optimization with ZOS-API.

Uses the thick-lens paraxial formula (in air) to enforce EFL at a specified
refractive index. RMS spot radius is measured from real OpticStudio reports.
This is an on-axis, single-wavelength demonstration, not a production design.
"""
import argparse
import csv
import math
import re
from pathlib import Path

from scipy.optimize import brentq, minimize_scalar
from build_lens import connect, DEFAULT_ROOT

RMS_RE = re.compile(r"RMS\s+Spot\s+Radius\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*(µm|μm|um|mm)", re.I)
# N-BK7 refractive index at 0.55 um, approximate catalog value. Confirm in
# OpticStudio glass catalog before reporting EFL as verified.
N_BK7_550 = 1.518522
T = 5.0
SEMI = 5.0
TARGET_EFL = 50.0


def read_report(path):
    raw = path.read_bytes()
    if raw.startswith((b'\xff\xfe', b'\xfe\xff')):
        return raw.decode('utf-16')
    if b'\x00' in raw[:200]:
        return raw.decode('utf-16-le')
    return raw.decode('utf-8-sig', errors='replace')


def rms_spot(system, zosapi, path):
    analysis = system.Analyses.New_Analysis(zosapi.Analysis.AnalysisIDM.StandardSpot)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.unlink(missing_ok=True)
        analysis.ApplyAndWaitForCompletion()
        analysis.GetResults().GetTextFile(str(path))
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError('Missing/empty Spot Diagram report: ' + str(path))
        match = RMS_RE.search(read_report(path))
        if not match:
            raise RuntimeError('Cannot parse RMS Spot Radius: ' + str(path))
        val = float(match.group(1))
        if match.group(2).lower() == 'mm':
            val *= 1000.0
        if not math.isfinite(val) or val <= 0:
            raise RuntimeError('Invalid RMS value: ' + str(val))
        return val
    finally:
        analysis.Close()


def power(r1, r2, n):
    # Thick lens in air; front and back surface powers, t = center thickness.
    p1 = (n - 1.0) / r1
    p2 = (1.0 - n) / r2
    return p1 + p2 - (T / n) * p1 * p2


def r2_for_efl(r1, efl, n):
    p1 = (n - 1.0) / r1
    # 1/efl = p1 + p2*(1 - T*p1/n)
    p2 = (1.0 / efl - p1) / (1.0 - T * p1 / n)
    if abs(p2) < 1e-12:
        raise ValueError('Rear surface is effectively plane')
    return (1.0 - n) / p2


def bfl_estimate(r1, r2, n):
    # Distance from rear vertex to paraxial focus for object at infinity.
    return (1.0 - T * ((n - 1.0) / r1) / n) / power(r1, r2, n)


def edge_thickness(r1, r2):
    # Restrict to biconvex positive-front, negative-rear radii.
    if r1 <= SEMI or r2 >= -SEMI or abs(r2) <= SEMI:
        return -math.inf
    s1 = r1 - math.sqrt(r1*r1 - SEMI*SEMI)
    s2 = abs(r2) - math.sqrt(r2*r2 - SEMI*SEMI)
    return T - s1 - s2


def symmetric_baseline_radius(efl, n):
    return brentq(lambda r: power(r, -r, n) - 1.0 / efl, 20.0, 200.0)


def build(system):
    system.New(False)
    lde = system.LDE
    lde.InsertNewSurfaceAt(2)
    obj, front, rear, image = (lde.GetSurfaceAt(i) for i in range(4))
    obj.Thickness = 1e10
    front.Comment = 'N-BK7 front / stop'
    front.Thickness = T
    front.Material = 'N-BK7'
    rear.Comment = 'Rear spherical surface'
    image.Comment = 'Image plane'
    system.SystemData.Aperture.ApertureValue = 10.0
    system.SystemData.Wavelengths.GetWavelength(1).Wavelength = 0.55
    field = system.SystemData.Fields.GetField(1)
    field.X = 0.0
    field.Y = 0.0
    return front, rear


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--zemax-root', default=DEFAULT_ROOT)
    p.add_argument('--outdir', default='results/v4')
    p.add_argument('--efl', type=float, default=TARGET_EFL)
    p.add_argument('--n', type=float, default=N_BK7_550,
                   help='N-BK7 index at 0.55 um; verify against glass catalog')
    p.add_argument('--min-edge', type=float, default=1.0)
    p.add_argument('--maxiter', type=int, default=28)
    args = p.parse_args()
    if args.efl <= 0 or args.n <= 1 or args.min_edge <= 0:
        p.error('Require positive EFL, n > 1 and positive minimum edge thickness')
    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    app = None
    rows = []
    best = {'rms': math.inf}
    count = 0

    def write_csv():
        with (out / 'spot_sweep.csv').open('w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['id','stage','r1_mm','r2_mm','bfl_mm','paraxial_efl_mm',
                        'edge_mm','rms_um','status','report'])
            w.writerows(rows)

    try:
        app, zosapi = connect(args.zemax_root)
        system = app.PrimarySystem
        front, rear = build(system)

        def measure(stage, r1, bfl):
            nonlocal count
            count += 1
            r2 = r2_for_efl(r1, args.efl, args.n)
            edge = edge_thickness(r1, r2)
            if edge < args.min_edge or bfl <= 0:
                return 1e8
            front.Radius = float(r1)
            rear.Radius = float(r2)
            rear.Thickness = float(bfl)
            path = out / 'reports' / f'{count:05d}_{stage}.txt'
            try:
                rms = rms_spot(system, zosapi, path)
                status = 'ok'
            except Exception as exc:
                rms = math.inf
                status = 'error: ' + str(exc)[:160]
            rows.append([count, stage, r1, r2, bfl, 1/power(r1,r2,args.n),
                         edge, rms, status, str(path.relative_to(out))])
            if math.isfinite(rms) and rms < best['rms']:
                best.update(rms=rms, r1=r1, r2=r2, bfl=bfl)
                system.SaveAs(str(out / 'best_so_far.zos'))
            if stage == 'baseline' or count % 15 == 0:
                print(f'{count:4d} {stage:9s} R1={r1:.3f} R2={r2:.3f} '
                      f'BFL={bfl:.3f} RMS={rms:.4f} um', flush=True)
            return rms if math.isfinite(rms) else 1e8

        r0 = symmetric_baseline_radius(args.efl, args.n)
        r20 = r2_for_efl(r0, args.efl, args.n)
        if edge_thickness(r0, r20) < args.min_edge:
            raise RuntimeError('Baseline violates minimum edge thickness')
        # Baseline has exactly the same EFL as all optimization candidates.
        bfl0 = bfl_estimate(r0, r20, args.n)
        baseline = measure('baseline', r0, bfl0)
        if baseline >= 1e8:
            raise RuntimeError('Baseline measurement failed; inspect reports')
        system.SaveAs(str(out / 'baseline.zos'))

        def optimize_focus(r1):
            try:
                r2 = r2_for_efl(r1, args.efl, args.n)
                if edge_thickness(r1, r2) < args.min_edge:
                    return 1e8
                guess = bfl_estimate(r1, r2, args.n)
                if not 5 < guess < 150:
                    return 1e8
                result = minimize_scalar(
                    lambda bfl: measure('focus', r1, bfl),
                    bounds=(max(1.0, guess-7.0), guess+7.0), method='bounded',
                    options={'xatol': 0.015, 'maxiter': 32})
                return float(result.fun)
            except (ValueError, ZeroDivisionError):
                return 1e8

        # Only positive-front / negative-rear biconvex shapes in this study.
        # Broad 35..95 mm interval; constraint is EFL, not independent R2.
        # Coarse multi-start scan protects against a poor 1-D local search.
        candidates = []
        for r1 in [35, 40, 45, 50, 55, 60, 70, 80, 95]:
            score = optimize_focus(float(r1))
            candidates.append((score, float(r1)))
        candidates.sort()
        for _, center in candidates[:3]:
            lo, hi = max(35.0, center-10), min(95.0, center+10)
            if hi-lo < 1:
                continue
            minimize_scalar(optimize_focus, bounds=(lo,hi), method='bounded',
                            options={'xatol': 0.1, 'maxiter': args.maxiter})
        if not math.isfinite(best['rms']):
            raise RuntimeError('No valid measured candidate')
        front.Radius = float(best['r1'])
        rear.Radius = float(best['r2'])
        rear.Thickness = float(best['bfl'])
        system.SaveAs(str(out / 'optimized_candidate.zos'))
        improvement = 100.0*(baseline-best['rms'])/baseline
        summary = (
            'V4 fixed-EFL, on-axis, 550 nm spherical N-BK7 singlet\n'
            f'Target paraxial EFL: {args.efl:.6f} mm\n'
            f'Assumed refractive index: {args.n:.8f}\n'
            'EPD: 10 mm; center thickness: 5 mm\n'
            f'Baseline: R1={r0:.6f} mm R2={r20:.6f} mm BFL={bfl0:.6f} mm\n'
            f'Baseline measured RMS: {baseline:.6f} um\n'
            f'Best measured RMS: {best["rms"]:.6f} um\n'
            f'Improvement: {improvement:.2f}%\n'
            f'Best R1: {best["r1"]:.6f} mm\n'
            f'Best R2: {best["r2"]:.6f} mm\n'
            f'Best BFL: {best["bfl"]:.6f} mm\n'
            f'Best edge thickness: {edge_thickness(best["r1"],best["r2"]):.6f} mm\n'
            'CAUTION: EFL enforced by analytical thick-lens formula using assumed n.\n'
            'VERIFY EFL against Zemax Prescription Data / cardinal points.\n'
            'MTF, Zernike, wavefront, tolerancing NOT automated or validated in V4.\n'
        )
        (out / 'summary.txt').write_text(summary, encoding='utf-8')
        print(summary)
        print('Saved results to', out)
    finally:
        write_csv()
        if app is not None:
            app.CloseApplication()


if __name__ == '__main__':
    main()
