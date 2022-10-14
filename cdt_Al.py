import os
from pylab import *

ds = [2.7, 0.027, 0.0027]
ts = [0.025, 10., 20., 30, 40., 48.0]

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'cd%dt%dAl'%(i,j)
        c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)            
        print(c)
        os.system(c)
        
