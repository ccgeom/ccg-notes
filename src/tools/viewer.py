import tools.converter as cvt
import pyvista as pv

mesh = pv.read('data/doubletorus.vtu')
mesh.plot(screenshot='data/doubletorus.png')
