import numpy as np
import pyvista as pv
import sys


def build_index(mesh, pionts):
    n_points = mesh.n_points
    neighbors = {}
    for i in range(n_points):
        print(i, end=':')
        sys.stdout.flush()
        if i not in neighbors.keys():
            neighbors[i] = []

        p = pionts[i]
        for j in range(n_points):
            if j not in neighbors.keys():
                neighbors[j] = []

            q = pionts[j]
            if i != j:
                d = np.sqrt(np.sum((p - q) * (p - q)))
                if d < 0.06:
                    neighbors[i].append((j, (q - p) / d))
                    neighbors[j].append((i, (p - q) / d))
                    print(len(neighbors[i]), end=',')
        print()

    return neighbors


def step(neighbors, pvalues):
    counter = 0
    result = np.copy(pvalues)
    for i in range(mesh.n_points):
        pval = pvalues[i]
        if pval > 0.8:                         # 如果 p 点正在燃烧
            if pval * 0.9 < 0.8:
                result[i] = 0.5
            else:
                result[i] = pval * 0.9

        vectors = []
        for j, v in neighbors[i]:
            qval = pvalues[j]
            if pval > 0.8:                      # 如果 p 点正在燃烧
                if qval == 0.7:                 # 如果 q 点尚未燃烧过
                    result[j] = 1.0             # 则点燃 q 点
                    counter += 1
            elif pval == 0.7:                   # 如果 p 点尚未燃烧过
                if qval > 0.9:                  # 如果 q 点刚刚燃烧
                    vectors.append(v)           # 则记录来火方向

        meeted = 0
        total = 0
        if len(vectors) > 0:
            for v1 in vectors:
                for v2 in vectors:
                    total += 1
                    if np.sum(v1 * v2) < 0:
                        meeted += 1
            if meeted * 3 > total:                       # 满足相遇条件
                result[i] = 0.0                 # 则 p 点是相遇点
            else:
                result[0] = 0.0                 # 强制把色阶拉回去的 workaround

    return counter, result


if __name__ == '__main__':
    print('reading mesh...')
    mesh = pv.read('data/doubletorus.vtu')
    mesh.point_arrays['pvalues'] = np.ones([mesh.n_points]) * 0.7  # 绿色的森林

    print('copy vertex...')
    pionts = np.array(mesh.points).copy()
    print('build index...')
    neighbors = build_index(mesh, pionts)

    print('initialize...')
    pvalues = np.ones([mesh.n_points]) * 0.7  # 绿色的森林
    pvalues[0] = 0.0                          # 强制把色阶拉回去的 workaround
    pvalues[5000] = 1.0                       # 初始着火处

    print('evolving...')
    for ix in range(50):
        mesh.point_arrays['pvalues'][:] = pvalues
        mesh.plot(scalars='pvalues', screenshot='data/firing_%02d.png' % ix, interactive=False)
        counter, result = step(neighbors, pvalues)
        pvalues = result
        print(counter)
