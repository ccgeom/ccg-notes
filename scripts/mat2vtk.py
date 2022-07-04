from glob import glob
from tools.converter import mat2vtk

if __name__ == '__main__':
    for mfile in glob('data/*.mat'):
        vfile = mfile.replace('.mat', '.vtk')
        mat2vtk(mfile, vfile)
        print(mfile, '->', vfile)
    print('finished')
