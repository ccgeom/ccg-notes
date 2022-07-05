import pyvista as pv

mesh = pv.read('data/eight.vtk')
mesh.plot(screenshot='data/eight.png')
