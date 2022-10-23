import os
from pylab import *
from multiprocessing import Pool

nproc = 16
ds = [0.234, 2.34]
ts = [0.025, 5.0, 10., 15., 20., 25., 30., 100., 300., 1e3, 3e3, 10e3]
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/dtB/d%dt%02dB'%(i,j)
    c = 'python dedx.py --zt=5 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)

p = Pool(processes=nproc)
p.map(run1dt, range(nd*nt))
    
