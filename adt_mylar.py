import os, sys
from pylab import *

ds = [1.35, 0.00864, 0.00432]
ts = [0.025, 3.0, 10.0]
i0 = 0
i1 = len(ds)-1
if len(sys.argv) > 1:
    i0 = int(sys.argv[1])
    i1 = i0
    if len(sys.argv) > 2:
        i1 = int(sys.argv[2])

for i in range(i0,i1+1):
    od = 'tmp/dt%dMylar'%(i)
    c = "python dedx.py --zc='1,6,8' --wc='4,5,2' --aa=2 --mloss=24 --d=%g --t=%g --od=%s"%(ds[i], ts[i], od)            
    print(c)
    os.system(c)
        
