from glob import glob
from tools.converter import m2mat

if __name__ == '__main__':
    for mfile in glob('data/*.m'):
        matfile = mfile.replace('.m', '.mat')
        m2mat(mfile, matfile)
        print(mfile, '->', matfile)
    print('finished')
