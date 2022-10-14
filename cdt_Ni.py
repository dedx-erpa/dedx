import os
from pylab import *

ds = [8.91, 0.0891, 0.00891]
ts = [0.025, 10., 20., 30., 42.0]

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'cd%dt%dNi'%(i,j)
        c = 'python dedx.py --zt=28 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
        print(c)
        os.system(c)
        
