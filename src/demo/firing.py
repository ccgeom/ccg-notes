import numpy as np
import pyvista as pv


def build_index(mesh, pionts):
    neighbors = {}
    for i in range(mesh.n_points):
        if i not in neighbors.keys():
            neighbors[i] = []

        p = pionts[i]
        for j in range(mesh.n_points):
            if j not in neighbors.keys():
                neighbors[j] = []

            q = pionts[j]
            if i != j:
                d = np.sqrt(np.sum((p - q) * (p - q)))
                if d < 0.05:
                    neighbors[i] = (j, (q - p) / d)
                    neighbors[j] = (i, (p - q) / d)

    return neighbors


def step(neighbors, pvalues):
    counter = 0
    for i in range(mesh.n_points):
        pval = pvalues[i]
        vectors = []
        for j, v in neighbors[i]:
            qval = pvalues[j]
            if pval > 0.75:                      # 如果 p 点正在燃烧
                if qval == 0.7:
                    pvalues[j] = 1.0             # 则点燃 q 点
                    counter += 1
            elif pval == 0.7:                    # 如果 p 点尚未燃烧过
                if qval > 0.75:                  # 如果 q 点正在燃烧
                    vectors.append(v)            # 则记录来火方向
        meet_cond = False
        if len(vectors) > 0:
            for v in vectors:
                if np.sum(v * vectors[0]) < 0:
                    meet_cond = True
                    break
            if meet_cond:                        # 满足相遇条件
                pvalues[i] = 0.0                 # 则 p 点是相遇点
            else:
                pvalues[0] = 0.0                 # 强制把色阶拉回去的 workaround

    return counter


if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')
    mesh.point_arrays['pvalues'] = np.ones([mesh.n_points]) * 0.7  # 绿色的森林

    print('build index...')
    pionts = np.array(mesh.points).copy()
    neighbors = build_index(mesh, pionts)

    print('initialize...')
    pvalues = np.ones([mesh.n_points]) * 0.7  # 绿色的森林
    pvalues[0] = 0.0                          # 强制把色阶拉回去的 workaround
    pvalues[5000] = 1.0                       # 初始着火处

    print('evolving...')
    for ix in range(13):
        mesh.point_arrays['pvalues'][:] = pvalues
        mesh.plot(scalars='pvalues', screenshot='data/firing_%02d.png' % ix, interactive=False)
        print(step())
