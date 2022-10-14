import os
from pylab import *

ds = [0.0891, 0.891, 8.91]
ts = 10**linspace(-1,4.,35)

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'd%dt%02dNi'%(i,j)
        c = 'python dedx.py --zt=28 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
            
        print(c)
        os.system(c)
        
