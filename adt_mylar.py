import os
from pylab import *

ds = [1.35, 0.00864, 0.00432]
ts = [0.025, 3.0, 10.0]

for i in range(1,len(ds)):
    od = 'tmp/dt%dMylar'%(i)
    c = "python dedx.py --zc='1,6,8' --wc='4,5,2' --aa=2 --mloss=24 --d=%g --t=%g --od=%s"%(ds[i], ts[i], od)            
    print(c)
    os.system(c)
        
