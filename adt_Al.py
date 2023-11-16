import os, sys
from pylab import *

ds = [2.7, 0.01728, 0.00864]
ts = [0.025, 4.5, 15.0]
i0 = 0
i1 = len(ds)-1
if len(sys.argv) > 1:
    i0 = int(sys.argv[1])
    i1 = i0
    if len(sys.argv) > 2:
        i1 = int(sys.argv[2])

for i in range(i0,i1+1):
    od = 'tmp/dt%dAl'%(i)
    c = 'python dedx.py --zt=13 --aa=2 --mloss=24 --d=%g --t=%g --od=%s'%(ds[i], ts[i], od)            
    print(c)
    os.system(c)
        
