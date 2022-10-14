import os
import numpy
from multiprocessing import Pool
from pfac import fac
import rd
from dts import *

def run_dedx(z, d, t, aa=2, zp=1,
             mloss=24, floss='NONE', floss0='NONE', bqp=-1E12):
    if type(z) == type(0):
        a = fac.ATOMICSYMBOL[z]
    else:
        a = z
        zc,wc = rd.zwc(z)
        if len(zc) == 1:
            z = zc[0]
            a = fac.ATOMICSYMBOL[z]
        else:
            z = 0
    odir = 'data/%s'%(a)
    if z > 0:
        c = 'python dedx.py --zt=%d --d=%g --t=%g --aa=%d --bqp=%-7.1E --mloss=%d --floss=%s --floss0=%s --od=%s --zp=%g'%(z, d, t, aa, bqp, mloss, floss, floss0, odir, zp)
    else:
        c = 'python dedx.py --fc=%s --d=%g --t=%g --aa=%d --bqp=%-7.1E --mloss=%d --floss=%s --floss0=%s --od=%s --zp=%g'%(a, d, t, aa, bqp, mloss, floss, floss0, odir, zp)
        
    print(c)
    os.system(c)

def run1zp(z, m):
    d = getden(z)
    if d > 0:
        t = tmin   
        run_dedx(z, d, t, aa=2, bqp=-1E12, floss='tct.dat', mloss=24)

def run1z(z):
    run1zp(z, -1)
    
def run_loop(xs):
    for a in xs:
        run1zp(a[0], a[1])
        
def dist_zs(izs, np):
    n = len(izs)

    xs = []
    for i in range(len(izs)):
        xs.append((izs[i],-1))
                
    a = [None]*np
    for i in range(0, n, np):
        ni = min(i+np, n)
        for j in range(i,ni):
            if a[j-i] is None:
                a[j-i] = []
            a[j-i].append(xs[j])
    return a

def run_azs(izs, np=16):
    a = dist_zs(izs, np)
    p = Pool(processes=np)
    p.map(run_loop, a)
    
    
