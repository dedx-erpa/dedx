import os
from pylab import *
from multiprocessing import Pool

nproc = 16
ds = array([0.01, 0.1, 1.0])*19.311
ts = [10**(-1.0+i*5.0/34) for i in range(35)]
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/dtAu/d%dt%02dAu'%(i,j)
    c = 'python dedx.py --zt=79 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)

if __name__ == '__main__':
    p = Pool(processes=nproc)
    p.map(run1dt, range(nd*nt))

