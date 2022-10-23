import os
from pylab import *

ds = [2.7, 0.01728, 0.00864]
ts = [0.025, 4.5, 15.0]

for i in range(len(ds)):
    od = 'tmp/dt%dAl'%(i)
    c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[i], od)            
    print(c)
    os.system(c)
        
