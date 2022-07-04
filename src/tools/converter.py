import numpy as np
import scipy.io
import vedo


def mread(mfile):
    vertexes, faces = [], []
    remap = {}
    with open(mfile) as f:
        for lnn, lnt in enumerate(f):
            fields = lnt.strip().replace('\t', ' ').replace('  ', ' ').split(' ')
            if len(fields) == 0:
                continue
            elif fields[0] == '#':
                continue
            else:
                typ, idx, d1, d2, d3 = fields
                if typ == 'Vertex':
                    vertexes.append([float(d1), float(d2), float(d3)])
                    remap[int(idx)] = len(vertexes)
                elif typ == 'Face':
                    faces.append([remap[int(d1)], remap[int(d2)], remap[int(d3)]])
    return vertexes, faces


def m2mat(mfile, matfile):
    vertexes, faces = mread(mfile)
    scipy.io.savemat(matfile, {'node': vertexes, 'elem': faces})


def m2vtk(mfile, vfile):
    vertexes, faces = mread(mfile)
    node = np.array(vertexes, dtype=np.float64)
    cell = np.array(faces, dtype=np.int_) - 1
    mesh = vedo.Mesh([node, cell])
    vedo.io.write(mesh, vfile)


def mat2vtk(matfile, vfile):
    data = scipy.io.loadmat(matfile)
    node = np.array(data['node'], dtype=np.float64)
    cell = np.array(data['elem'] - 1, dtype=np.int_)
    mesh = vedo.Mesh([node, cell])
    vedo.io.write(mesh, vfile)




