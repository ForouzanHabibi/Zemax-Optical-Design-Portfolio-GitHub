"""Build a reproducible sequential singlet in OpticStudio 2024 R1 via ZOS-API.
No existing lens files are modified. Windows + pythonnet + licensed OpticStudio required.
"""
from pathlib import Path
import argparse
import sys

DEFAULT_ROOT = r'C:\Program Files\Ansys Zemax OpticStudio 2024 R1.00'

def connect(root):
    if sys.platform != 'win32':
        raise RuntimeError('This script must run on Windows with OpticStudio installed.')
    import clr
    folder = Path(root)
    helper = folder / 'ZOSAPI_NetHelper.dll'
    if not helper.is_file():
        raise FileNotFoundError(helper)
    clr.AddReference(str(helper))
    import ZOSAPI_NetHelper
    if not ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize():
        raise RuntimeError('ZOS-API initialization failed')
    install = Path(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())
    clr.AddReference(str(install / 'ZOSAPI.dll'))
    clr.AddReference(str(install / 'ZOSAPI_Interfaces.dll'))
    import ZOSAPI
    connection = ZOSAPI.ZOSAPI_Connection()
    app = connection.CreateNewApplication()
    if app is None or not app.IsValidLicenseForAPI:
        raise RuntimeError('No valid ZOS-API license')
    return app, ZOSAPI

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--zemax-root', default=DEFAULT_ROOT)
    p.add_argument('--output', default='results/singlet_baseline.zos')
    args = p.parse_args()
    app = None
    try:
        app, zosapi = connect(args.zemax_root)
        system = app.PrimarySystem
        system.New(False)
        # Sequential: OBJECT, STOP, two refracting surfaces, IMAGE.
        lde = system.LDE
        lde.InsertNewSurfaceAt(2)
        lde.InsertNewSurfaceAt(3)
        obj = lde.GetSurfaceAt(0)
        front = lde.GetSurfaceAt(1)
        back = lde.GetSurfaceAt(2)
        image = lde.GetSurfaceAt(3)
        obj.Thickness = 1.0e10  # object at practical infinity (mm)
        front.Comment = 'Spherical singlet front / STOP'
        front.Radius = 50.0
        front.Thickness = 5.0
        front.Material = 'N-BK7'
        back.Comment = 'Spherical singlet rear'
        back.Radius = -50.0
        back.Thickness = 95.0  # initial image distance; optimize later
        image.Comment = 'Image plane'
        # Wavelength and field are intentionally left at OpticStudio defaults.
        # Set them explicitly and validate in the GUI before claiming performance.
        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        system.SaveAs(str(output))
        print('Saved:', output)
        print('Inspect in OpticStudio: system units, wavelengths, aperture, fields, stop, spot and MTF.')
        print('This is an initial lens model, NOT an optimized or validated design.')
    finally:
        if app is not None:
            app.CloseApplication()

if __name__ == '__main__':
    main()
