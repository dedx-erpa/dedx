import os
from pylab import *

ds = [0.234, 2.34]
ts = [0.025, 5.0, 10., 15., 20., 25., 30., 100., 300., 1e3, 3e3, 10e3]

for i in range(len(ds)):
    for j in range(len(ts)):
        od = 'd%dt%02dB'%(i,j)
        if j == 0:
            c = 'python dedx.py --zt=5 --aa=2 --mloss=24 --floss=tst2.dat --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
        else:
            c = 'python dedx.py --zt=5 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[j], od)
            
        print(c)
        os.system(c)
        
