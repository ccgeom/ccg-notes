
import sys

import numpy as np
import pyvista as pv


class Manifold:
    def __init__(self):
        pass


class Firefront:
    def __init__(self, cutlocus):
        pass

    def cut(self, point):
        pass

    def join(self, point):
        pass

    def leave(self, point):
        pass

    def step(self):
        pass


class FireBulk:
    def __init__(self):
        pass

    def add(self, point):
        pass

    def minus(self, point):
        pass

    def step(self):
        pass


class CutLocus:
    def __init__(self):
        pass

    def join(self, point):
        pass


class WildFireMethod:
    def __init__(self, manifold):
        self.manifold = manifold
        self.cutlocus = CutLocus()
        self.firefront = Firefront(self.cutlocus)
        self.firebulk = FireBulk()

    def step(self):
        for face in self.firefront:
            for fngbr in self.manifold.neighbor_faces(face):
                if not self.firebulk.contains(fngbr):
                    self.firefront.add(fngbr)
            self.firefront.minus(face)
            self.firebulk.add(face)


def make_dual(points, cells, n_cells):
    ix, cur = 0, 0
    centroid = []
    cells_begin = []
    cells_end = []
    while ix < n_cells:
        sz = cells[cur]
        assert sz > 2
        ps = points[cells[cur + 1:cur + sz]]
        assert ps.shape[1] == sz
        centroid.append(np.mean(ps, axis=1))
        cells_begin.append(cur + 1)
        cells_end.append(cur + sz)
        cur = cur + sz + 1
        ix += 1
    assert cur == cells.shape[0]
    return np.array(centroid), np.array(cells_begin), np.array(cells_end)


def build_index(points, n_points):
    neighbors = {}
    for i in range(n_points):
        print(i, end=':')
        sys.stdout.flush()
        if i not in neighbors.keys():
            neighbors[i] = set()

        p = points[i]
        for j in range(n_points):
            if j not in neighbors.keys():
                neighbors[j] = set()

            q = points[j]
            if i != j:
                d = np.sqrt(np.sum((p - q) * (p - q)))
                if d < 0.01:
                    neighbors[i].add((j, (q - p) / d))
                    neighbors[j].add((i, (p - q) / d))
                    print(len(neighbors[i]), end=',')
        print()

    return neighbors


def step(points, dual, cells, cells_begin, cells_end, pneighbors, cneighbors, cvalues, last_fireline, last_cutlines):
    counter = 0
    fireline = []
    n_points, n_cells = points.shape[0], dual.shape[0]
    result = np.copy(cvalues)
    for i in range(n_cells):
        ppos = dual[i]
        pval = cvalues[i]
        if pval > 0.8:                         # 如果 p 点正在燃烧
            if pval * 0.9 < 0.8:               # 如果 p 点马上烧尽
                result[i] = 0.5
            else:
                result[i] = pval * 0.9         # 如果 p 点未烧尽

        if pval == 0.7:                        # 如果 p 点尚未燃烧过
            for j, v in cneighbors[i]:
                qpos = dual[j]
                qval = cvalues[j]
                if pval > 0.9:                 # 如果 p 点正在燃烧
                    if qval == 0.7:            # 如果 q 点尚未燃烧过
                        result[j] = 1.0        # 则点燃 q 点
                        counter += 1

    return counter, result, fireline, cutlines



if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')
    n_points = mesh.n_points
    n_cells = mesh.n_cells
    mesh.cell_arrays['cvalues'] = np.zeros([n_cells])
    mesh.point_arrays['pvalues'] = np.zeros([n_points])

    print('copy data...')
    points = np.array(mesh.points).copy()
    cells = np.array(mesh.cells).copy()

    print('constructing dual...')
    dual, cells_begin, cells_end = make_dual(points, cells, mesh.n_cells)

    print('build index...')
    pneighbors = build_index(points, mesh.n_cells)
    cneighbors = build_index(dual, mesh.n_cells)

    print('initialize...')
    firefront = Firefront()
    cutline = Cutline()

    cvalues = np.ones([mesh.n_cells]) * 0.7   # 绿色的森林
    cvalues[0] = 0.0                          # 强制把色阶拉回去的 workaround
    cvalues[10000] = 1.0                      # 初始着火处
    path = None

    print('evolving...')
    counter, ix = 1, 0
    while counter != 0:
        print(ix, counter)
        mesh.cell_arrays['cvalues'][:] = cvalues

        plt = pv.Plotter()
        plt.add_mesh(mesh, scalars='cvalues')
        if path is not None:
            plt.add_mesh(mesh.extract_points(path), color='black')

        plt.show(screenshot='data/fire_%03d.png' % ix, interactive=False)

        counter, cvalues, path = step(points, dual, cells, cells_begin, cells_end, cneighbors, pneighbors, cvalues, path, mesh.n_cells)
        ix += 1
