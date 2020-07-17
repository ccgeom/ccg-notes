
import sys
import itertools
import random

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
            self.n_faces = mesh.n_cells

            print('copy data...')
            self.points = np.array(mesh.points).copy()
            self.faces = np.array(mesh.cells).copy()

            print('constructing dual...')
            self.dual, self.faces_begin, self.faces_end = self.make_dual(self.points, self.faces, self.n_faces)

            print('indexing point2faces...')
            self.point2faces = self.index_point2faces(self.faces, self.n_faces, self.n_points)

            print('build neighbor index...')
            self.neighbor_points = self.build_neighbors(self.points, self.n_points, self.points, self.n_points, 0.05)
            self.neighbor_faces = self.build_neighbors(self.dual, self.n_faces, self.dual, self.n_faces, 0.05)
            self.neighbor_point2face = self.build_neighbors(self.points, self.n_points, self.dual, self.n_faces, 0.05)

            self.sort_orientation()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["mesh"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.mesh = None

    def index_point2faces(self, faces, n_faces, n_points):
        ix, cur = 0, 0
        point2faces = [[] for _ in range(n_points)]
        while ix < n_faces:
            sz = faces[cur]
            assert sz > 2
            for p in faces[cur + 1:cur + sz]:
                point2faces[p].append(ix)
            cur = cur + sz + 1
            ix += 1
        assert cur == faces.shape[0]
        return point2faces

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
            centroid.append(np.mean(ps, axis=0))
            faces_begin.append(cur + 1)
            faces_end.append(cur + sz)
            cur = cur + sz + 1
            ix += 1
        assert cur == faces.shape[0]
        return np.array(centroid), np.array(faces_begin), np.array(faces_end)

    def build_neighbors(self, points1, n_points1, points2, n_points2, thrshold):
        neighbors = [list() for _ in range(n_points1)]
        flag = len(points1) != len(points2) or (points1 != points2).any()
        for i in range(n_points1):
            print(i, end=':')
            sys.stdout.flush()

            p = points1[i]
            for j in range(n_points2):
                q = points2[j]
                e = p - q
                d = np.sqrt(np.sum(e * e))
                if d < thrshold:
                    if flag:
                        neighbors[i].append(j)
                    else:
                        if i != j:
                            neighbors[i].append(j)

            print(len(neighbors[i]), end=',')
            print()

        return neighbors

    def sort_orientation(self):
        checker = EvolvingBoundary(self)
        counter = 0
        face = 0
        faces = {ix: -1 for ix in range(self.n_faces)}
        checked = set()
        candidate = set()
        while len(checked) != self.n_faces:
            try:
                if face not in checked:
                    checker.begin()
                    ps = self.edges(face)
                    checker.insert(ps)
                    checked.add(face)
                    for p in ps:
                        fs = manifold.point2faces[p]
                        for f in fs:
                            candidate.add(f)
                    while len(candidate) > 0:
                        face = candidate.pop()

                    checker.commit()
            except AssertionError:
                index = faces[face]
                begin, end = self.faces_begin[face], self.faces_end[face] + 1
                points = self.faces[begin:end]
                permute = list(itertools.permutations(points))
                self.faces[begin:end] = permute[index]
                faces[face] -= 1
                counter += 1
                assert faces[face] + len(permute) > 0
        print(checker.cycles)
        return counter

    def edges(self, face):
        begin = self.faces_begin[face]
        end = self.faces_end[face]
        return self.faces[begin:end+1]


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
                if a < b:
                    key, sgn = (a, b), 1
                if a > b:
                    key, sgn = (b, a), -1
                if key not in self.workarea:
                    self.workarea[key] = 0
                self.workarea[key] = self.workarea[key] + sgn

    def commit(self):
        next = {}
        for a, b in self.workarea.keys():
            val = self.workarea[(a, b)]
            if val > 0:
                assert val == 1
                next[a] = b
            if val < 0:
                assert val == -1
                next[b] = a
        paths = []
        while len(next) > 0:
            path = []
            key = next.keys().__iter__().__next__()
            while key in next:
                path.append(key)
                temp = next[key]
                next.pop(key)
                key = temp
            paths.append(np.array(path, dtype=np.int))
        self.cycles = paths
        self.workarea.clear()

    def insert(self, edges):
        sz = len(edges)
        assert sz > 2
        for ix in range(sz):
            a = edges[ix]
            b = edges[(ix + 1) % sz]
            if a < b:
                key, sgn = (a, b), 1
            if a > b:
                key, sgn = (b, a), -1
            if key not in self.workarea:
                self.workarea[key] = 0

            workarea = self.workarea.copy()
            workarea[key] = workarea[key] + sgn
            assert np.abs(workarea[key]) < 2.0
            self.workarea = workarea


class Firefront(EvolvingBoundary):
    def __init__(self, manifold, start_face):
        super().__init__(manifold)
        self.begin()
        self.insert(self.manifold.edges(start_face))
        self.commit()


class Firetail(EvolvingBoundary):
    def __init__(self, manifold):
        super().__init__(manifold)


class FireBulk:
    name = 'firebulk'
    def __init__(self, manifold, start_face):
        self.manifold = manifold
        self.bulk = set()
        self.listeners = []

        self.values = np.ones([manifold.mesh.n_cells]) * 0.7   # 绿色的森林
        self.values[0] = 0.0                                   # 强制把色阶拉回去的 workaround
        self.values[start_face] = 1.0                          # 初始着火处
        manifold.mesh.cell_arrays[self.name] = self.values

    def add_vanish_listener(self, lstnr):
        self.listeners.append(lstnr)

    def on_vanish(self, face):
        for lstnr in self.listeners:
            lstnr.on_firebulk_vanished(face)

    def add(self, face):
        self.values[face] = 1.0
        self.bulk.add(face)

    def minus(self, face):
        self.bulk.remove(face)
        self.on_vanish(face)

    def contains(self, face):
        return face in self.bulk

    def step(self):
        for face in self.bulk.copy():
            val = self.values[face]
            if val * 0.9 < 0.7:
                self.values[face] = 0.5
                self.minus(face)
            else:
                self.values[face] = val * 0.9


class CutLocus:
    def __init__(self, manifold):
        self.manifold = manifold
        self.points = set()
        self.path = []

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
        self.firebulk.add_vanish_listener(self)

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
        firefront_cycles = self.firefront.cycles
        self.begin()
        for cycle in firefront_cycles:
            for point in cycle:
                for fngbr in self.manifold.neighbor_point2face[point]:
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
        manifold = pickle.load(mf.open(mode='rb'))
        if manifold.sort_orientation() > 0:
            pickle.dump(manifold, mf.open(mode='wb'))
        mesh = pv.read('data/doubletorus.vtu')
        manifold.mesh = mesh
    else:
        print('reading mesh...')
        mesh = pv.read('data/doubletorus.vtu')
        print('constructing manifold...')
        manifold = Manifold(mesh)
        print('writing manifold...')
        pickle.dump(manifold, mf.open(mode='wb'))

    wildfire = WildFireSweepingMethod(manifold, start_face=10000)

    print('evolving...')
    counter, ix = 1, 1
    while counter != 0:
        print(ix, counter)
        counter = wildfire.step()

        plt = pv.Plotter()
        manifold.mesh.cell_arrays[wildfire.firebulk.name][:] = wildfire.firebulk.values
        plt.add_mesh(manifold.mesh, scalars=wildfire.firebulk.name)

        firefront_cycles = wildfire.firefront.cycles
        for cycle in firefront_cycles:
            plt.add_mesh(manifold.mesh.extract_points(cycle), color='red')
        cutlocus = wildfire.cutlocus.path
        if len(cutlocus) > 1:
            plt.add_mesh(manifold.mesh.extract_points(cutlocus), color='black')

        plt.show(screenshot='data/wildfire_%03d.png' % ix, interactive=False)
        ix += 1
