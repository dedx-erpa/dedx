import os
import numpy

ts = 10**numpy.linspace(-1.,4.,21)
for i in range(len(ts)):
    f = 't%02d.dat'%i
    if os.path.exists(f):
        continue
    c = 'python dief.py %s %7.2f 0'%(f, ts[i])
    print(c)
    os.system(c)
