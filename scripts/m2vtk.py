from glob import glob
from tools.converter import m2vtk

if __name__ == '__main__':
    for mfile in glob('data/*.m'):
        vfile = mfile.replace('.m', '.vtk')
        m2vtk(mfile, vfile)
        print(mfile, '->', vfile)
    print('finished')
