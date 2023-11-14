import os
from pylab import *
from multiprocessing import Pool

nproc = 16
ds = [0.0891, 0.891, 8.91]
ts = [10**(-1.0+i*5.0/34) for i in range(35)]
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/dtNi/d%dt%02dNi'%(i,j)
    c = 'python dedx.py --zt=28 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)

p = Pool(processes=nproc)
p.map(run1dt, range(nd*nt))

