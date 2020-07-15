
import sys

import numpy as np
import pyvista as pv


# 在下同调的语言里，引入一个差分算子 \delta = x_i - x_{i-1}，用来刻画演化问题。
# 在火烧法里，引入记号 f 表示火头线， B 表示火体面，t 表示火尾线.
# f_i = \delta \partial \delta^+ B_i, f_0 = 0
# t_i = \delta \partial \delta^- B_i, t_0 = 0
# \delta^+ B_i 表示进入火体的面
# \delta^- B_i 表示离开火体的面

class Manifold:
    def __init__(self, mesh=None):
        if mesh != None:
            self.mesh = mesh
            self.n_points = mesh.n_points
            self.n_faces = mesh.n_faces

            print('copy data...')
            self.points = np.array(mesh.points).copy()
            self.faces = np.array(mesh.faces).copy()

            print('constructing dual...')
            self.dual, self.faces_begin, self.faces_end = self.make_dual(self.points, self.faces, self.n_faces)

            print('build neighbor index...')
            self.neighbor_points = self.build_neighbors(self.points, self.n_points, 0.05)
            self.neighbor_faces = self.build_neighbors(self.dual, self.n_faces, 0.05)

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

    def make_dual(self, points, faces, n_faces):
        ix, cur = 0, 0
        centroid = []
        faces_begin = []
        faces_end = []
        while ix < n_faces:
            sz = faces[cur]
            assert sz > 2
            ps = points[faces[cur + 1:cur + sz]]
            assert ps.shape[1] == sz
            centroid.append(np.mean(ps, axis=1))
            faces_begin.append(cur + 1)
            faces_end.append(cur + sz)
            cur = cur + sz + 1
            ix += 1
        assert cur == faces.shape[0]
        return np.array(centroid), np.array(faces_begin), np.array(faces_end)

    def edges(self, face):
        return self.faces[self.faces_begin[face]:self.faces_end[face]]


class EvolvingBoundary:
    def __init__(self, manifold):
        self.manifold = manifold
        self.cycles = []
        self.workarea = {}

    def begin(self):
        self.workarea.clear()
        for c in self.cycles:
            r = c[::-1]
            l = len(r)
            for ix in range(l):
                a = r[ix - 1]
                b = r[ix]
                c = r[(ix + 1) % l]
                if a not in self.workarea:
                    self.workarea[a] = 0, b
                if b not in self.workarea:
                    self.workarea[b] = 0, c

                self.workarea[a][0] = self.workarea[a][0] - 1
                self.workarea[b][0] = self.workarea[b][0] + 1

    def commit(self):
        pass

    def insert(self, edges):
        for ix in len(edges):
            a = edges[ix - 1]
            if a not in self.workarea:
                self.workarea[a] = 0
            b = edges[ix]
            if b not in self.workarea:
                self.workarea[b] = 0

            self.workarea[a] = self.workarea[a] - 1
            self.workarea[b] = self.workarea[b] + 1


class Firefront(EvolvingBoundary):
    def __init__(self, manifold, start_face):
        super(manifold)
        self.listeners = []

    def add_meet_listener(self, lstnr):
        self.listeners.append(lstnr)

    def on_meet(self, point):
        for lstnr in self.listeners:
            lstnr.on_firefront_meeted(self, point)


class Firetail(EvolvingBoundary):
    def __init__(self, manifold):
        super(manifold)


class FireBulk:
    name = 'firebulk'
    def __init__(self, manifold, start_face):
        self.manifold = manifold
        self.bulk = set()
        self.listeners = []

        self.values = np.ones([manifold.mesh.n_faces]) * 0.7   # 绿色的森林
        self.values[0] = 0.0                                   # 强制把色阶拉回去的 workaround
        self.values[start_face] = 1.0                          # 初始着火处
        manifold.mesh.face_arrays[self.name] = self.values

    def add_vanish_listener(self, lstnr):
        self.listeners.append(lstnr)

    def on_vanish(self, face):
        for lstnr in self.listeners:
            lstnr.on_firetail_vanished(self, face)

    def add(self, face):
        self.values[face] = 1.0
        self.bulk.add(face)

    def minus(self, face):
        self.values[face] = 0.5
        self.bulk.remove(face)
        self.on_vanish(face)

    def contains(self, face):
        return face in self.bulk

    def step(self):
        for face in self.bulk:
            val = self.values[face]
            if val * 0.9 < 0.7:
                self.minus(face)
            else:
                self.values[face] = val * 0.9


class CutLocus:
    def __init__(self, manifold):
        self.manifold = manifold
        self.points = set()
        self.path = None

    def add(self, point):
        self.points.add(point)
        self.path = list(self.points)

    def minus(self, point):
        self.points.remove(point)
        self.path = list(self.points)


class WildFireSweepingMethod:
    def __init__(self, manifold, start_face):
        self.manifold = manifold
        self.start_face = start_face
        self.cutlocus = CutLocus(manifold)
        self.firebulk = FireBulk(manifold, start_face)
        self.firefront = Firefront(manifold, start_face)
        self.firetail = Firetail(manifold)
        self.firefront.add_meet_listener(self)
        self.firebulk.add_vanish_listener(self)

    def on_firefront_meeted(self, point):
        pass

    def on_firebulk_vanished(self, face):
        self.firetail.insert(manifold.edges(face))

    def begin(self):
        self.firefront.begin()
        self.firetail.begin()

    def commit(self):
        self.firefront.commit()
        self.firetail.commit()

    def step(self):
        counter = 0
        firefront_path = self.firefront.path
        self.begin()
        for point in firefront_path:
            for fngbr in self.manifold.neighbor_faces_by_point(point):
                if not self.firebulk.contains(fngbr):
                    self.firebulk.add(fngbr)
                    self.firefront.insert(manifold.edges(fngbr))
                    counter += 1
        self.commit()
        self.firebulk.step()
        return counter


if __name__ == '__main__':
    import pathlib
    import pickle

    mf = pathlib.Path('data/doubletorus.pkl')
    if mf.exists():
        print('reading manifold...')
        manifold = pickle.load(mf)
    else:
        print('reading mesh...')
        mesh = pv.read('data/doubletorus.vtu')
        print('constructing manifold...')
        manifold = Manifold(mesh)
        print('writing manifold...')
        pickle.dump(manifold, mf.open(mode='w'))

    wildfire = WildFireSweepingMethod(manifold, start_face=5000)

    print('evolving...')
    counter, ix = 1, 0
    while counter != 0:
        print(ix, counter)
        counter = wildfire.step()

        plt = pv.Plotter()
        plt.add_mesh(manifold.mesh, scalars=wildfire.firebulk.name)

        firefront = wildfire.firefront.path
        if firefront is not None:
            plt.add_mesh(manifold.mesh.extract_points(firefront), color='red')
        cutlocus = wildfire.cutlocus.path
        if cutlocus is not None:
            plt.add_mesh(manifold.mesh.extract_points(cutlocus), color='black')

        plt.show(screenshot='data/wildfire_%03d.png' % ix, interactive=False)
        ix += 1
