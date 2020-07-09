import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

data = scipy.io.loadmat('data/doubletorus.mat')
X = data['X']
Y = data['Y']
U = data['U']

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, U, rstride=1, cstride=1, cmap=cm.jet)
fig.colorbar(surf)

plt.show()