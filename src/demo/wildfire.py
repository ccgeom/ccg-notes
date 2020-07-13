
import sys

import numpy as np
import pyvista as pv


class Manifold:
    def __init__(self, mesh):
        self.mesh = mesh
        self.n_points = mesh.n_points
        self.n_cells = mesh.n_cells
        mesh.point_arrays['pvalues'] = np.zeros([self.n_points])

        print('copy data...')
        self.points = np.array(mesh.points).copy()
        self.cells = np.array(mesh.cells).copy()

        print('constructing dual...')
        self.dual, self.cells_begin, self.cells_end = self.make_dual(self.points, self.cells, self.n_cells)

        print('build index...')
        self.neighbor_points = self.build_neighbors(self.points, self.n_points)
        self.neighbor_faces = self.build_neighbors(self.dual, self.n_cells)

    def build_neighbors(points, n_points):
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


class Firefront:
    def __init__(self, manifold):
        self.manifold = manifold
        self.listeners = []
        self.path = None

    def add_meet_listener(self, lstnr):
        self.listeners.append(lstnr)

    def on_meet(self, point):
        for lstnr in self.listeners:
            lstnr.on_firefront_meeted(self, point)


class FireBulk:
    def __init__(self, manifold):
        self.manifold = manifold
        self.values = np.ones([manifold.mesh.n_cells]) * 0.7   # 绿色的森林
        self.values[0] = 0.0                                   # 强制把色阶拉回去的 workaround
        self.values[10000] = 1.0                               # 初始着火处
        manifold.mesh.cell_arrays['firebulk'] = self.values
        self.bulk = set()

    def add(self, cell):
        self.values[cell] = 1.0
        self.bulk.add(cell)

    def minus(self, cell):
        self.values[cell] = 0.5
        self.bulk.remove(cell)

    def step(self):
        for cell in self.bulk:
            val = self.values[cell]
            if val * 0.9 < 0.7:
                self.minus(cell)
            else:
                self.values[cell] = val * 0.9


class CutLocus:
    def __init__(self, manifold):
        self.manifold = manifold

    def join(self, point):
        pass

    def on_firefront_meeted(self, point):
        pass


class WildFireMethod:
    def __init__(self, manifold):
        self.manifold = manifold
        self.cutlocus = CutLocus(manifold)
        self.firefront = Firefront(manifold)
        self.firebulk = FireBulk(manifold)
        self.firefront.add_meet_listener(self.cutlocus)

    def step(self):
        counter = 0
        for face in self.firefront:
            for fngbr in self.manifold.neighbor_faces(face):
                if not self.firebulk.contains(fngbr):
                    self.firefront.add(fngbr)
                    counter += 1
            self.firefront.minus(face)
            self.firebulk.add(face)
        return counter


if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')

    print('initialize...')
    manifold = Manifold(mesh)
    wildfire = WildFireMethod(manifold)

    print('evolving...')
    counter, ix = 1, 0
    while counter != 0:
        print(ix, counter)
        counter = wildfire.step()

        plt = pv.Plotter()
        plt.add_mesh(mesh, scalars='firebulk')

        firefront = wildfire.firefront
        if firefront.path is not None:
            plt.add_mesh(mesh.extract_points(firefront.path), color='red')
        cutlocus = wildfire.cutlocus
        if cutlocus.path is not None:
            plt.add_mesh(mesh.extract_points(cutlocus.path), color='black')

        plt.show(screenshot='data/wildfire_%03d.png' % ix, interactive=False)
        ix += 1
