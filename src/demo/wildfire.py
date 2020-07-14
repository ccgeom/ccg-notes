
import sys

import numpy as np
import pyvista as pv


class Manifold:
    def __init__(self, mesh):
        self.mesh = mesh
        self.n_points = mesh.n_points
        self.n_cells = mesh.n_cells

        print('copy data...')
        self.points = np.array(mesh.points).copy()
        self.cells = np.array(mesh.cells).copy()

        print('constructing dual...')
        self.dual, self.cells_begin, self.cells_end = self.make_dual(self.points, self.cells, self.n_cells)

        print('build neighbor index...')
        #self.neighbor_points = self.build_neighbors(self.points, self.n_points, 0.05)
        #self.neighbor_faces = self.build_neighbors(self.dual, self.n_cells, 0.05)

    def build_neighbors(self, points, n_points, thrshold):
        neighbors = []
        for i in range(n_points):
            print(i, end=':')
            sys.stdout.flush()
            if i >= len(neighbors):
                neighbors.append(set())

            p = points[i]
            for j in range(n_points):
                if j >= len(neighbors):
                    neighbors.append(set())

                q = points[j]
                if i != j:
                    d = np.sqrt(np.sum((p - q) * (p - q)))
                    if d < thrshold:
                        neighbors[i].add(j)
                        neighbors[j].add(i)
                        print(len(neighbors[i]), end=',')
            print()

        return neighbors

    def make_dual(self, points, cells, n_cells):
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


class Firetail:
    def __init__(self, manifold):
        self.manifold = manifold
        self.listeners = []
        self.path = None

    def add_vanish_listener(self, lstnr):
        self.listeners.append(lstnr)

    def on_vanish(self, point):
        for lstnr in self.listeners:
            lstnr.on_firetail_vanished(self, point)


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
                self.values[cell] = 0.5
            else:
                self.values[cell] = val * 0.9


class CutLocus:
    def __init__(self, manifold):
        self.manifold = manifold

    def join(self, point):
        pass


class WildFireSweepingMethod:
    def __init__(self, manifold):
        self.manifold = manifold
        self.cutlocus = CutLocus(manifold)
        self.firefront = Firefront(manifold)
        self.firetail = Firetail(manifold)
        self.firebulk = FireBulk(manifold)
        self.firefront.add_meet_listener(self)
        self.firetail.add_vanish_listener(self)

    def on_firefront_meeted(self, point):
        pass

    def on_firetail_vanished(self, point):
        pass

    def step(self):
        counter = 0
        for face in self.firebulk.faces:

        for point in self.firefront.path:
            for fngbr in self.manifold.neighbor_faces_by_point(point):
                if not self.firebulk.contains(fngbr):
                    self.firebulk.add(fngbr)
                    counter += 1
        return counter


if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')

    print('initialize...')
    manifold = Manifold(mesh)
    wildfire = WildFireSweepingMethod(manifold)

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
