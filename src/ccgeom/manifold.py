
import numpy as np

# following the idea of halfedge data structure on David Gu's lecture
# https://www3.cs.stonybrook.edu/~gu/lectures/lecture_8_halfedge_data_structure.pdf
# and adapt it to a numpy friendly representatives
# * vertexes: all position of each vertex
# * faces: all position of each centroid of faces
# * edges: all position of each centroid of edges
# * halfedges: all vectors of each halfedge
# * vertexes2vertexes


class Manifold:
    def __init__(self, vtk_mesh=None):
        if vtk_mesh != None:
            self.mesh = vtk_mesh # a VTK mesh structure

            self.n_vertexes = vtk_mesh.n_points
            self.n_faces = vtk_mesh.n_cells

            cells = np.array(self.mesh.cells).copy()
            self.vertexes = np.array(self.mesh.points).copy()
            self.faces, cells_begin, cells_end = make_dual(self.n_faces, self.vertexes, cells)
            self.edges, self.halfedges = make_edges(self.n_faces, self.vertexes, cells, cells_begin, cells_end)
            self.n_edges = self.edges.shape[0]
            self.n_halfedges = self.halfedges.shape[0]

            self.adjacency_vertexes = None
            self.adjacency_faces = None
            self.adjacency_edges = None
            self.adjacency_halfedges = None

            self.adjacency_vertexes2faces = None
            self.adjacency_vertexes2edges = None
            self.adjacency_vertexes2halfedges = None

            self.adjacency_faces2vertexes = None
            self.adjacency_faces2edges = None
            self.adjacency_faces2halfedges = None

            self.adjacency_edges2vertexes = None
            self.adjacency_edges2faces = None
            self.adjacency_edges2halfedges = None

            self.adjacency_halfedges2vertexes = None
            self.adjacency_halfedges2faces = None
            self.adjacency_halfedges2edges = None

            self.orientation = None

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["mesh"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.mesh = None


def make_dual(n_faces, points, faces):
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


def make_edges(n_faces, vertexes, cells, cells_begin, cells_end):
    total = 0
    for ix in range(n_faces):
        begin, end = cells_begin[ix], cells_end[ix]
        sz = end - begin + 1
        total += sz

    cur = 0
    edges = np.array([total, 3], dtype=np.float64)
    halfedges = np.array([2 * total, 3], dtype=np.float64)
    for ix in range(n_faces):
        begin, end = cells_begin[ix], cells_end[ix]
        sz = end - begin + 1
        pxs = vertexes[cells[begin:end]]
        for p in range(sz):
            src, tgt = pxs[p - 1], pxs[p]
            edges[cur] = (src + tgt) / 2
            halfedges[2 * cur] = tgt - src
            halfedges[2 * cur + 1] = src - tgt
            cur += 1

    return edges, halfedges