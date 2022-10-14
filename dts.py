import numpy, os
from pfac import fac

fa = fac.ATOMICSYMBOL
ds = list(numpy.loadtxt('dens.txt', unpack=1, usecols=1))
ss = list(numpy.loadtxt('dens.txt', unpack=1, usecols=0, dtype=str))
zs = []
for i in range(len(ss)):
    if ss[i].isdigit():
        zs.append(int(ss[i]))
        ss[i] = fa[zs[-1]]
    else:
        zs.append(0)

ts = numpy.arange(0.0, 1000.1, 50.)
tmin = 0.025
ts[0] = tmin
mts = 10**numpy.linspace(-1.,4.,21)

def getden(x):
    if type(x) == type(0):
        return d4z(x)
    else:
        return d4a(x)
    
def d4z(z):    
    try:
        iz = zs.index(z)
    except:
        print('cannot find z for rho: %d'%z)
        return 0.0
    return ds[iz]

def d4a(a):
    try:
        iz = ss.index(a)
    except:
        print('cannot find a for rho: %s'%a)
        return 0.0
    return ds[iz]
