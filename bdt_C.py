import os
from pylab import *
from multiprocessing import Pool

nproc = 12
ds = [0.05, 0.1, 0.5]
ts = [0.025, 10.,  20., 30.]

nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/d%dt%02dC'%(i,j)
    c = 'python dedx.py --zt=6 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
    
    print(c)
    os.system(c)
        

p = Pool(processes=nproc)
p.map(run1dt, range(nd*nt))
