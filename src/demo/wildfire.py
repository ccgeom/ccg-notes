
import sys

import numpy as np
import pyvista as pv


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


def build_index(pionts, n_points):
    neighbors = {}
    for i in range(n_points):
        print(i, end=':')
        sys.stdout.flush()
        if i not in neighbors.keys():
            neighbors[i] = set()

        p = pionts[i]
        for j in range(n_points):
            if j not in neighbors.keys():
                neighbors[j] = set()

            q = pionts[j]
            if i != j:
                d = np.sqrt(np.sum((p - q) * (p - q)))
                if d < 0.01:
                    neighbors[i].add(j)
                    neighbors[j].add(i)
                    print(len(neighbors[i]), end=',')
        print()

    return neighbors


def step(cells, cells_begin, cells_end, cneighbors, cvalues, path, n_cells):
    counter = 0
    result = np.copy(cvalues)
    for i in range(n_cells):
        pval = cvalues[i]
        if pval > 0.8:                         # 如果 p 点正在燃烧
            if pval * 0.9 < 0.8:               # 如果 p 点马上烧尽
                result[i] = 0.5
            else:
                result[i] = pval * 0.9         # 如果 p 点未烧尽

        for j in cneighbors[i]:
            qval = cvalues[j]
            if pval > 0.8:                     # 如果 p 点正在燃烧
                if qval == 0.7:                # 如果 q 点尚未燃烧过
                    result[j] = 1.0            # 则点燃 q 点
                    counter += 1
                elif qval > 0.8:
                    ppoints = cells[cells_begin[i]:cells_end[i]]
                    qpoints = cells[cells_begin[j]:cells_end[j]]
                    commons = set(ppoints).intersection(set(qpoints))
                    if path is None:
                        path = np.array(list(commons))
                    else:
                        path = np.concatenate([path, np.array(list(commons))])

    return counter, result, path


if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')
    mesh.cell_arrays['cvalues'] = np.zeros([mesh.n_cells])
    mesh.point_arrays['pvalues'] = np.zeros([mesh.n_points])

    print('copy data...')
    pionts = np.array(mesh.points).copy()
    cells = np.array(mesh.cells).copy()

    print('constructing dual...')
    dual, cells_begin, cells_end = make_dual(pionts, cells, mesh.n_cells)

    print('build index...')
    neighbors = build_index(dual, mesh.n_cells)

    print('initialize...')
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
        plt.show(screenshot='data/fire_%03d.png' % ix, interactive=False)

        counter, result, path = step(cells, cells_begin, cells_end, neighbors, cvalues, path, mesh.n_cells)
        if path is not None:
            plt.add_mesh(mesh.extract_points(path), color='red')
        cvalues = result
        ix += 1
