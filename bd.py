import os
from pylab import *
from multiprocessing import Pool

nproc = 6
ds = [4e0, 4e1, 4e2, 4e3, 4e4, 4e5]
ts = [2e3]
nd = len(ds)
nt = len(ts)

def run1dt(ij):
    i = int(ij/nt)
    j = ij%nt
    od = 'tmp/d%dB'%i
    c = 'python dedx.py --zt=5 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
    print(c)
    os.system(c)

if __name__ == '__main__':
    p = Pool(processes=nproc)
    p.map(run1dt, range(nd*nt))
    
