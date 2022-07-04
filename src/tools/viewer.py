from vedo import *

mesh = Mesh('data/eight.vtk')
mesh.c('white').lighting('glossy')

p1 = Point([0, +1, 0], c='yellow')
p2 = Point([0, -1, 0], c='yellow')

l1 = Light(p1, c='yellow')
l2 = Light(p2, c='yellow')

show(mesh, l1, l2, p1, p2, "eight", axes=True)
io.screenshot('eight.png')
