import os
import numpy
from multiprocessing import Pool
from pfac import fac
from pfac import aa as aam
import rd
from dts import *

def run_dedx(z, d, t, aa=2, zp=1, npaa=0, maa=1,
             mloss=24, floss='NONE', floss0='NONE', bqp=-1E12):
    if type(z) == type(0):
        a = fac.ATOMICSYMBOL[z]
    else:
        a = z
        zc,wc = aam.zw4c(z)
        if len(zc) == 1:
            z = zc[0]
            a = fac.ATOMICSYMBOL[z]
        else:
            z = 0

    odir = 'data/%s'%a
    if mloss == 0:
        frho = odir
        odir = '%s0'%odir
        
    if z > 0:
        c = 'python dedx.py --zt=%d --d=%g --t=%g --aa=%d --bqp=%-7.1E --mloss=%d --floss=%s --floss0=%s --od=%s --zp=%g'%(z, d, t, aa, bqp, mloss, floss, floss0, odir, zp)
    else:
        c = 'python dedx.py --fc=%s --d=%g --t=%g --aa=%d --bqp=%-7.1E --mloss=%d --floss=%s --floss0=%s --od=%s --zp=%g --npaa=%d --maa=%d'%(a, d, t, aa, bqp, mloss, floss, floss0, odir, zp, npaa, maa)

    if mloss == 0:
        c = '%s --frho=%s/rho.functions'%(c,frho)
        
    print(c)
    os.system(c)

def run1zp(z, npaa=0, maa=1, aa=2, mloss=24):
    d = getden(z)
    if d > 0:
        t = tmin   
        run_dedx(z, d, t, floss='tct.dat', npaa=npaa, maa=maa,
                 aa=aa, mloss=mloss)
    
def run_loop(xs):
    if xs is None:
        return
    for a in xs:
        run1zp(a[0], aa=a[1], mloss=a[2], npaa=1, maa=1)
        
def dist_zs(izs, np, aa=2, mloss=24):
    n = len(izs)

    xs = []
    for i in range(len(izs)):
        xs.append((izs[i],aa,mloss))
                
    a = [None]*np
    for i in range(0, n, np):
        ni = min(i+np, n)
        for j in range(i,ni):
            if a[j-i] is None:
                a[j-i] = []
            a[j-i].append(xs[j])
    return a

def run_azs(izs, np=16, aa=2, mloss=24):
    a = dist_zs(izs, np, aa=aa, mloss=24)
    p = Pool(processes=np)
    p.map(run_loop, a)
    
    
