import random
import sys
import itertools

import numpy as np
import pyvista as pv

from scipy.spatial.distance import cdist
from alive_progress import alive_bar


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
            self.center = mesh.center
            self.area = mesh.area
            self.volume = mesh.volume
            self.bounds = mesh.bounds

            print('copy points...')
            self.points = np.array(mesh.points).copy()

            print('constructing faces...')
            self.faces = np.array(mesh.faces).copy()
            self.face_normals = np.array(mesh.face_normals).copy()
            self.faces_centroid, self.faces_begin, self.faces_end = self.build_faces(self.points, self.faces, self.n_faces)

            print('indexing point2faces...')
            self.point2faces = self.index_point2faces(self.faces, self.n_faces, self.n_points)

            print('build neighbor index...')
            with alive_bar(3) as bar:
                self.neighbor_points = self.build_neighbors(self.points, self.n_points, self.points, self.n_points, 0.05)
                bar()
                self.neighbor_faces = self.build_neighbors(self.faces_centroid, self.n_faces, self.faces_centroid, self.n_faces, 0.05)
                bar()
                self.neighbor_point2face = self.build_neighbors(self.points, self.n_points, self.faces_centroid, self.n_faces, 0.05)
                bar()

            #print('sort orientation...')
            #self.sort_orientation()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["mesh"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.mesh = None

    def edges(self, face):
        begin = self.faces_begin[face]
        end = self.faces_end[face]
        return self.faces[begin:end+1]

    def index_point2faces(self, faces, n_faces, n_points):
        ix, cur = 0, 0
        point2faces = [[] for _ in range(n_points)]

        with alive_bar(n_faces) as bar:
            while ix < n_faces:
                sz = faces[cur]
                assert sz > 2
                for p in faces[cur + 1:cur + sz]:
                    point2faces[p].append(ix)
                cur = cur + sz + 1
                ix += 1
                bar()
        assert cur == faces.shape[0]
        return point2faces

    def build_faces(self, points, faces, n_faces):
        ix, cur = 0, 0
        centroid = []
        faces_begin = []
        faces_end = []
        with alive_bar(n_faces) as bar:
            while ix < n_faces:
                sz = faces[cur]
                assert sz > 2
                ps = points[faces[cur + 1:cur + sz]]
                assert ps.shape[1] == sz
                centroid.append(np.mean(ps, axis=0))
                faces_begin.append(cur + 1)
                faces_end.append(cur + sz)
                assert np.max(faces[cur + 1:cur + sz]) < self.n_points
                cur = cur + sz + 1
                ix += 1
                bar()
        assert cur == faces.shape[0]
        return np.array(centroid, dtype=np.float32), np.array(faces_begin, dtype=np.int64), np.array(faces_end, dtype=np.int64)

    def build_neighbors(self, points1, n_points1, points2, n_points2, thrshold):
        dmatrix = cdist(points1, points2, metric='euclidean')
        neighbors = [[] for _ in range(n_points1)]
        flag = len(points1) != len(points2) or (points1 != points2).any()
        for i in range(n_points1):
            for j in range(n_points2):
                d = dmatrix[i, j]
                if d < thrshold:
                    if flag:
                        neighbors[i].append(j)
                    else:
                        if i != j:
                            neighbors[i].append(j)
        return neighbors

    def sort_orientation(self):
        checker = EvolvingBoundary(self)
        checked = set()
        candidate = set([0])
        flips = np.zeros(self.n_faces, dtype=np.int8)
        with alive_bar(self.n_faces) as bar:
            while len(checked) < self.n_faces:
                while len(candidate) > 0:
                    face = candidate.pop()
                    try:
                        ps = self.edges(face)
                        checker.begin()
                        checker.insert(ps)
                        checker.commit()
                        checked.add(face)
                        bar()
                    except AssertionError:
                        begin, end = self.faces_begin[face], self.faces_end[face] + 1
                        reversed(self.faces[begin:end])
                        flips[face] += 1
                    for p in ps:
                        fs = self.point2faces[p]
                        for f in fs:
                            if f not in checked:
                                candidate.add(f)

        print(checker.cycles)
        return np.sum(flips)


class EvolvingBoundary:
    def __init__(self, manifold):
        self.manifold = manifold
        self.cycles = []
        self.working_area = {}

    def begin(self):
        self.working_area.clear()
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
                if key not in self.working_area:
                    self.working_area[key] = 0
                self.working_area[key] = self.working_area[key] + sgn

    def commit(self):
        next = {}
        for (a, b), val in self.working_area.items():
            if val > 0:
                #assert val == 1
                next[a] = b
            if val < 0:
                #assert val == -1
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
            paths.append(np.array(path, dtype=np.int64))
        self.cycles = paths
        self.working_area.clear()

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

            if key not in self.working_area:
                self.working_area[key] = 0
            val = self.working_area[key] + sgn
            #assert np.abs(val) < 2.0
            self.working_area[key] = val


class FireFront(EvolvingBoundary):
    def __init__(self, manifold, start_face):
        super().__init__(manifold)
        self.begin()
        self.insert(self.manifold.edges(start_face))
        self.commit()


class FireTail(EvolvingBoundary):
    def __init__(self, manifold):
        super().__init__(manifold)


class FireBulk:
    name = 'firebulk'

    def __init__(self, manifold, start_face):
        self.manifold = manifold
        self.bulk = set()
        self.listeners = []

        self.values = np.ones([manifold.mesh.n_cells]) * 0.7   # 黄色的森林
        self.values[0] = 0.0                                   # 强制把色阶拉回去的 workaround
        self.values[start_face] = 1.0                          # 初始着火处
        manifold.mesh.cell_data[self.name] = self.values

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
        self.values[0] = 0.0   # 强制把色阶拉回去的 workaround


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
        self.firefront = FireFront(manifold, start_face)
        self.firetail = FireTail(manifold)
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

    vtk_file = 'data/eight.vtk'
    pkl_file = 'data/eight.pkl'

    mf = pathlib.Path(pkl_file)
    if mf.exists():
        print('reading manifold...')
        manifold = pickle.load(mf.open(mode='rb'))
        #print('check orientation...')
        #if manifold.sort_orientation() > 0:
        #    pickle.dump(manifold, mf.open(mode='wb'))
        mesh = pv.read(vtk_file)
        manifold.mesh = mesh
    else:
        print('reading mesh...')
        mesh = pv.read(vtk_file)
        print('constructing manifold...')
        manifold = Manifold(mesh)
        print('writing manifold...')
        pickle.dump(manifold, mf.open(mode='wb'))

    wildfire = WildFireSweepingMethod(manifold, start_face=3000)

    print('evolving...')
    counter, ix = 1, 1
    while counter != 0:
        print(ix, counter)
        counter = wildfire.step()

        plt = pv.Plotter()
        manifold.mesh.cell_data[wildfire.firebulk.name][:] = wildfire.firebulk.values
        plt.add_mesh(manifold.mesh, scalars=wildfire.firebulk.name)

        firefront_cycles = wildfire.firefront.cycles
        for cycle in firefront_cycles:
            plt.add_mesh(manifold.mesh.extract_points(cycle), color='red')
        cutlocus = wildfire.cutlocus.path
        if len(cutlocus) > 1:
            plt.add_mesh(manifold.mesh.extract_points(cutlocus), color='black')

        plt.show(screenshot='data/wildfire_%03d.png' % ix, interactive=False)
        ix += 1
