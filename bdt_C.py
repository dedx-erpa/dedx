import os
from pylab import *

ds = [0.5]
ts = [10.,  20., 30.]

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'd%dt%02dC'%(i,j)
        c = 'python dedx.py --zt=6 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
            
        print(c)
        os.system(c)
        
