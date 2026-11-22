"""V5: fixed-curvature, fixed-EFL aspheric A4 + image-plane optimization.

Windows, OpticStudio 2024 R1, ZOS-API. Uses actual Zemax Zernike RMS reports.
No claims about MTF, tolerances or manufacture are made automatically.
"""
import argparse
import csv
import math
import re
from pathlib import Path
from scipy.optimize import minimize
from build_lens import connect, DEFAULT_ROOT
from optimize_singlet_v4 import build, rms_spot, edge_thickness, N_BK7_550, T, SEMI, power

RMS_WAVE = re.compile(r'From integration of the rays\s*:\s*RMS\s*\(to chief\)\s*:\s*([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s*waves', re.I)


def read_text(path):
    data = path.read_bytes()
    if data.startswith((b'\xff\xfe', b'\xfe\xff')):
        return data.decode('utf-16')
    if b'\x00' in data[:200]:
        return data.decode('utf-16-le')
    return data.decode('utf-8-sig', errors='replace')


def zernike_rms(system, zosapi, report):
    names = ('ZernikeStandardCoefficients', 'ZernikeStandard')
    enum = zosapi.Analysis.AnalysisIDM
    kind = next((getattr(enum, n) for n in names if hasattr(enum, n)), None)
    if kind is None:
        raise RuntimeError('Zernike Standard analysis enum not found in this ZOS-API version')
    analysis = system.Analyses.New_Analysis(kind)
    try:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.unlink(missing_ok=True)
        analysis.ApplyAndWaitForCompletion()
        analysis.GetResults().GetTextFile(str(report))
        if not report.exists() or not report.stat().st_size:
            raise RuntimeError('Zernike report was not generated')
        txt = read_text(report)
        m = RMS_WAVE.search(txt)
        if not m:
            raise RuntimeError('Cannot parse Zernike RMS; inspect ' + str(report))
        return float(m.group(1))
    finally:
        analysis.Close()


def asphere_sag(r, radius, a4):
    # Even asphere: conic constant 0; r in mm, A4 in mm^-3.
    c = 1.0 / radius
    rad = 1.0 - (c*r)**2
    if rad <= 0:
        return math.nan
    return c*r*r/(1 + math.sqrt(rad)) + a4*r**4


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--zemax-root', default=DEFAULT_ROOT)
    p.add_argument('--outdir', default='results/v5')
    p.add_argument('--r1', type=float, default=35.0)
    p.add_argument('--r2', type=float, default=-95.124393)
    p.add_argument('--bfl', type=float, default=47.150517)
    p.add_argument('--a4-bound', type=float, default=0.0001, help='A4 absolute bound in mm^-3')
    p.add_argument('--min-edge', type=float, default=1.0)
    p.add_argument('--maxiter', type=int, default=65)
    args = p.parse_args()
    if args.a4_bound <= 0 or args.bfl <= 0:
        p.error('Positive A4 bound and image distance required')
    out = Path(args.outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    app = None
    rows = []
    best = {'rms': math.inf}
    try:
        app, zosapi = connect(args.zemax_root)
        system = app.PrimarySystem
        front, rear = build(system)
        front.Radius = args.r1
        rear.Radius = args.r2
        rear.Thickness = args.bfl
        if edge_thickness(args.r1, args.r2) < args.min_edge:
            raise RuntimeError('Spherical baseline violates minimum edge thickness')
        system.SaveAs(str(out / 'spherical_reference.zos'))
        # Switch the FRONT surface from Standard to Even Asphere.
        lde = system.LDE
        surface_type = zosapi.Editors.LDE.SurfaceType.EvenAspheric
        settings = front.GetSurfaceTypeSettings(surface_type)
        front.ChangeType(settings)
        # Even Asphere: Parameter 1 = conic; Parameter 2 = 4th-order coefficient.
        # Verify these labels in the OpticStudio LDE for your installed version.
        col = zosapi.Editors.LDE.SurfaceColumn
        conic_cell = front.GetSurfaceCell(col.Par1)
        a4_cell = front.GetSurfaceCell(col.Par2)
        conic_cell.DoubleValue = 0.0
        a4_cell.DoubleValue = 0.0
        if abs(a4_cell.DoubleValue) > 1e-15:
            raise RuntimeError('Unable to initialize A4 parameter')
        system.SaveAs(str(out / 'asphere_zero_A4.zos'))

        def measure(stage, a4, bfl):
            a4, bfl = float(a4), float(bfl)
            if abs(a4) > args.a4_bound or not (args.bfl - 6 <= bfl <= args.bfl + 6):
                return 1e3
            # Front sag at edge minus rear sag contribution, using same sign convention.
            edge = T - asphere_sag(SEMI, args.r1, a4) + asphere_sag(SEMI, args.r2, 0.0)
            if not math.isfinite(edge) or edge < args.min_edge:
                return 1e3
            a4_cell.DoubleValue = a4
            rear.Thickness = bfl
            index = len(rows)
            path = out / 'reports' / f'{index:04d}_{stage}.txt'
            try:
                rms = zernike_rms(system, zosapi, path)
                if not math.isfinite(rms) or rms <= 0:
                    raise RuntimeError('Invalid RMS wavefront')
                status = 'ok'
            except Exception as e:
                rms = math.inf
                status = str(e)
            rows.append([stage, a4, bfl, edge, rms, status, str(path.relative_to(out))])
            if rms < best['rms']:
                best.update(rms=rms, a4=a4, bfl=bfl, edge=edge)
                system.SaveAs(str(out / 'best_so_far.zos'))
            if index % 10 == 0:
                print(f'{index:3d} {stage:10s} A4={a4:.8g} BFL={bfl:.5f} RMS={rms:.6f} waves', flush=True)
            return rms if math.isfinite(rms) else 1e3

        baseline = measure('spherical', 0, args.bfl)
        if baseline >= 1e3:
            raise RuntimeError('Zernike RMS baseline failed. Check first report in results/v5/reports')
        # Optimize a scaled A4 (units of 1e-5 mm^-3) and image distance.
        # Multistart guards against local minima; bounded by physical feasibility.
        starts = [(0.0, args.bfl), (-2.0, args.bfl), (2.0, args.bfl)]
        for a, f in starts:
            minimize(lambda x: measure('search', x[0]*1e-5, x[1]), [a, f],
                     method='Nelder-Mead',
                     options={'maxiter': args.maxiter, 'xatol': 1e-4, 'fatol': 1e-5})
        a4_cell.DoubleValue = best['a4']
        rear.Thickness = best['bfl']
        system.SaveAs(str(out / 'optimized_asphere.zos'))
        # Spot result is measured separately at the final wavefront-optimized focus.
        spot = rms_spot(system, zosapi, out / 'final_spot.txt')
        with (out / 'summary.txt').open('w', encoding='utf-8') as f:
            f.write('V5: even-asphere A4 and image-plane wavefront optimization\n')
            f.write(f'Fixed R1={args.r1:.6f} mm, R2={args.r2:.6f} mm; N-BK7, 550 nm, EPD=10 mm\n')
            f.write(f'Analytical paraxial EFL (vertex curvature): {1/power(args.r1,args.r2,N_BK7_550):.6f} mm\n')
            f.write(f'Baseline Zernike RMS at V4 spot focus: {baseline:.8f} waves\n')
            f.write(f'Best Zernike RMS: {best["rms"]:.8f} waves\n')
            f.write(f'Best A4: {best["a4"]:.10g} mm^-3\n')
            f.write(f'Best image distance: {best["bfl"]:.6f} mm\n')
            f.write(f'Best edge thickness: {best["edge"]:.6f} mm\n')
            f.write(f'Final measured RMS spot: {spot:.6f} um\n')
            f.write('VALIDATE: Even Asphere Par2=A4, EFL, Zernike Z11, MTF, sag/slope, tolerances in OpticStudio.\n')
            f.write('NOT VALIDATED: manufacturing feasibility and full-field performance.\n')
        print((out / 'summary.txt').read_text(encoding='utf-8'))
    finally:
        with (out / 'wavefront_sweep.csv').open('w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['stage', 'A4_mm_minus3', 'image_distance_mm', 'edge_mm', 'rms_wavefront_waves', 'status', 'report'])
            w.writerows(rows)
        if app is not None:
            app.CloseApplication()


if __name__ == '__main__':
    main()
