
import numpy as np
import scipy.io

from fealpy.mesh import TriangleMesh
from fealpy.writer import MeshWriter


def mat2vtk(mfile, vfile):
    data = scipy.io.loadmat(mfile)
    node = np.array(data['node'], dtype=np.float64)
    cell = np.array(data['elem'] - 1, dtype=np.int_)
    mesh = TriangleMesh(node, cell)
    writer = MeshWriter(mesh)
    writer.write(fname=vfile)




