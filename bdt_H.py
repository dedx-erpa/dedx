import os
from pylab import *

ds = [0.01]
ts = [1.]

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'd%dt%02dH'%(i,j)
        c = 'python dedx.py --zt=1 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
            
        print(c)
        os.system(c)
        
