
import numpy as np

# following the idea of halfedge data structure on David Gu's lecture
# https://www3.cs.stonybrook.edu/~gu/lectures/lecture_8_halfedge_data_structure.pdf
# and adapt it to a numpy friendly representatives
# * vertexes: all position of each vertex
# * faces: all position of each centroid of faces
# * edges: all position of each centroid of edges
# * halfedges: all vectors of each halfedge

class Manifold:
    def __init__(self, vtk_mesh=None):
        if vtk_mesh != None:
            self.mesh = vtk_mesh # a VTK mesh structure

            self.n_vertexes = vtk_mesh.n_points
            self.n_faces = vtk_mesh.n_cells
            self.n_edges = 0
            self.n_halfedges = 0

            self.vertexes = np.array(self.mesh.points).copy()

            cells = np.array(self.mesh.cells).copy()
            self.faces, cells_begin, cells_end = make_dual(self.n_faces, self.vertexes, cells)

            self.edges, self.halfedges = make_edges(self.n_faces, self.vertexes, cells)

            print('copy data...')
            self.vertexes = np.array(self.mesh.points).copy()
            self.faces = np.array(self.mesh.cells).copy()

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


def make_dual(points, faces, n_faces):
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


def make_halfedge():
    pass