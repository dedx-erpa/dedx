import os
from pylab import *
from multiprocessing import Pool

nproc = 16

ds = [2.7, 0.027, 0.0027]
ts = [0.025, 10., 20., 30, 40., 48.0]
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/cdtAl/d%dt%02dAl'%(i,j)
    c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)

p = Pool(processes=nproc)
p.map(run1dt, range(nd*nt))
        
