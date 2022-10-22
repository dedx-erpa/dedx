import os
from pylab import *
from multiprocessing import Pool

nproc = 16
ds = [0.027,0.27, 2.7]
ts = 10**linspace(-1,4.0,35)
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/dtAl/d%dt%02dAl'%(i,j)
    c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)
        

p = Pool(processes=nproc)
p.map(run1dt, range(nd*nt))
