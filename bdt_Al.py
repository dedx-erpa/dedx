import os
from pylab import *

ds = [0.027,0.27, 2.7]
ts = 10**linspace(-1,4.0,35)

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'd%dt%02dAl'%(i,j)
        c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
            
        print(c)
        os.system(c)
        
