
import numpy as np
import scipy.io

from fealpy.mesh import TriangleMesh
from fealpy.mesh import HalfEdgeMesh2d
from fealpy.geometry import SphereSurface
from fealpy.writer import MeshWriter



def mat2vtk(mfile):
    data = scipy.io.loadmat('data/doubletorus.mat')
    node = np.array(data['node'], dtype=np.float64)
    cell = np.array(data['elem'] - 1, dtype=np.int_)
    mesh = TriangleMesh(node, cell)


