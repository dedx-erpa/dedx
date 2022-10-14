from scipy import integrate, special, interpolate, optimize
from numpy import *
from pylab import *
from pfac import util, fac, consts
from multiprocessing import Pool
from time import time
import sys

def lastzero(y):
    w = where(y[1:]*y[:-1] <= 0)[0]
    if len(w) > 0:
        return w[-1]
    return len(y)

def exp1p(x):
    if x < 1:
        return 1/(1+exp(x))
    x = exp(-x)
    return x/(1+x)

def log1p(x):
    if (x > 10):
        return x+exp(-x)
    if (x < -10):
        return exp(x)
    return log(1+exp(x))

def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]

def cdiv(a, b):
    c = b[0]*b[0] + b[1]*b[1]
    d = cmul(a,(b[0],-b[1]))
    return d[0]/c,d[1]/c

def lfc(a, g):
    b = cmul(a, g)
    return cdiv(a, (1.0-b[0], -b[1]))

def elf(a):
    return a[1]/((1+a[0])**2+a[1]**2)

def wgint(y,x):
    return exp((y*y-1)*x*x);

def wgfun(x):
    if (x > 50 or x < 0.02):
        return 2*sqrt(pi)*x/(2*x**2+1)
    return integrate.quad(wgint, 0.0, 1.0, args=(x,))[0]*2*sqrt(pi)*x
        
def gint(y,te,ge,x):
    if (y == 0 or x == 0):
        return 0.0
    xyp = x+y
    xym = x-y
    t = (y*y - ge)/te
    yi = (y*exp1p(t))*log(abs((xyp)/(xym)))
    return yi

def fdint(t,x,n):
    return (t**n)*exp1p(t-x)

def subgint(x,te,ge,x0,x1):
    a = (te,ge,x)
    return integrate.quad(gint, x0, x1, args=a)[0]

def gfun(x, te, ge, eps=1e-10):
    xt2 = max(0.1*te,25*te+ge)
    xt = sqrt(xt2)
    tes = sqrt(abs(te))
    if x < xt:
        xm = 0.9*x
        xp = 1.1*x
        r0 = subgint(x, te, ge, 0.0, xm)
        r1 = subgint(x, te, ge, xm, x-eps)
        r2 = subgint(x, te, ge, x+eps, xp)
        r3 = subgint(x, te, ge, xp, xt)
        r4 = special.erfcx(xt/tes)
        if r4 > 0:
            r4 = exp(log(r4)+(ge-xt2)/te)
            r4 *= sqrt(pi)*tes*x
        r = r0+r1+r2+r3+r4
    else:
        r0 = subgint(x, te, ge, 0.0, xt)
        r4 = special.erfcx(xt/tes)*sqrt(pi)/2.0 + xt/tes
        if (r4 > 0):
            r4 = exp(log(r4)+(ge-xt2)/te)
            r4 *= tes*te/x
        r = r0+r4
    return r

def fdfun(x,n):
    x1 = max(0.1,25.0+x)
    r = integrate.quad(fdint, 0.0, x1, args=(x,n))[0]
    r1 = special.gamma(1+n)*special.gammaincc(1+n,x1)
    if (r1 > 0):
        r += exp(x+log(r1))
    return r

def interpgfun(x0, x1, nx, te, ge, rd=0):
    tes = sqrt(te)
    a = fdfun(ge/te,-0.5)*tes
    b = 2./3
    xa = sqrt(b/a)
    xs = 10**linspace(x0, x1, nx)*xa
    gs = array([gfun(t, te, ge) for t in xs])
    w = where(gs > 0)
    igf = interpolate.interp1d(log(xs[w]), log(gs[w]),
                               bounds_error=False,
                               fill_value = 'extrapolate')
    if rd == 0:
        return igf
    else:
        return igf, xs, gs, xa

def spar(de,te):
    de = de*(consts.RBOHR*1e-8)**3
    te = te/consts.HARTREE_EV
    rs = (4/3.*pi*de)**(-1./3)
    a = (9*pi/4.)**(-1./3)
    kf = 1.0/(a*rs)
    ef = 0.5*kf**2
    te = te/ef
    ge = util.FM1MP(4, (2./3)*te**(-1.5))*te
    wp = sqrt(4*pi*de)
    ga = 1/(te*ef*rs)
    return kf,de,te,ge,wp,rs,ga

def gapcorr(z, u, iu, kf):
    q = z*2
    d2 = (iu*iu - u*u)*2*q
    d = sqrt(d2)
    q2 = q*q
    q3 = q2*q
    q5 = q3*q2
    t0 = 1/q2
    t1 = -(arctan((2*q+q2)/d) + arctan((2*q-q2)/d))*(d/(2*q3))
    t2 = log((d2+(2*q+q2)**2)/(d2+(2*q-q2)**2))
    t2 *= (d2/(8*q5) + 1/(2*q3)-1/(8*q))
    a = (t0+t1+t2)*2/(kf*pi)            
    return a, 0.0
    
def rpa(z, u, iu, kf, te, ge, igfun):
    if (u < iu):
        return gapcorr(z, u, iu, kf)

    u = sqrt(u*u - iu*iu)
    uzp = u+z
    uzm = u-z
    uzma = max(1e-31,abs(uzm))
    if not igfun is None:
        igp = exp(igfun(log(uzp)))
        igm = exp(igfun(log(uzma)))
        if (uzm < 0):
            igm = -igm
    else:
        igp = gfun(uzp, te, ge)
        igm = gfun(uzma, te, ge)
        if (uzm < 0):
            igm = -igm
    z3kf = z**3*kf
    a = (1./(4*z3kf*pi))*(igp-igm)
    xp = (ge-uzp**2)/te
    xm = (ge-uzm**2)/te
    b = (te/(8*z3kf))*(log1p(xm)-log1p(xp))

    return a,b

def wpm(z0, u0, iu0, kf, te, ge, igf):
    if (u0 < iu0):
        return gapcorr(z0, u0, iu0, kf)
    if (iu0 > 0):
        u0 = sqrt(u0*u0 - iu0*iu0)
    te = sqrt(te*te + 0.4)
    tes = sqrt(te)
    u = u0/tes
    z = z0/tes
    x2 = 4/(3*sqrt(pi)*pi*kf*te*te)
    c1 = x2/(8*z**3)
    c2 = c1*pi
    uzp = u+z
    uzm = u-z
    uzma = abs(uzm)
    if not igf is None:
        igp = exp(igf(log(uzp)))
        igm = exp(igf(log(uzma)))
        if (uzm < 0):
            igm = -igm
    else:
        igp = wgfun(uzp)
        igm = wgfun(uzma)
        if (uzm < 0):
            igm = -igm
    a = c1*(igp - igm)
    b = c2*(exp(-uzm**2) - exp(-uzp**2))
    return a,b

def mermin(k,w,nu,r1,r0):
    a = cmul((w,nu),r1)
    b = cdiv(r1,r0)
    b = (w-nu*b[1],nu*b[0])
    return cdiv(a, b)
    
def testdief(de, te, rk, nu, nw=500, ig=0, iwpm=0):
    kf, de, te, ge, wp, rs, ga = spar(de, te)
    k = rk*kf
    wk = k*kf
    w = linspace(0.0, 4.0, nw)*wk
    u = w/(k*kf)
    iu = nu/(k*kf)

    if ig > 0:
        if (iwpm == 0):
            igf = interpgfun(-3., 3., 500, te, ge)
        else:
            igf = interpwgfun(-3., 3., 500, te, ge)
    else:
        igf = None
    if (iwpm == 0):
        r1 = array([rpa(0.5*rk, t, iu, kf, te, ge, igf) for t in u])
    else:
        r1 = array([wpm(0.5*rk, t, iu, kf, te, ge, igf) for t in u])
        
    r1 = (r1[:,0],r1[:,1])

    return kf,te,ge,wp,w,r1,elf(r1)

def giu(k, kf):
    x = k/kf
    A = 0.029
    a = (4/(9*pi))**(1/3.)
    rs = 1./(a*kf)
    z = 4*sqrt(a*rs/pi)
    b0 = 0.0621814
    b1 = 9.81379
    b2 = 2.82224
    b3 = 0.736411
    rs2 = rs*rs
    rs3 = rs2*rs
    rss = sqrt(rs)
    rs15 = rss*rs
    rs25 = rs2*rss
    rn = 1/rs3 + b1/rs25
    rd = 1 + b1*rss + b2*rs + b3*rs15
    rnd = -3.0/(rs3*rs) - 2.5*b1/(rs25*rs)
    rdd = 0.5*b1/rss + b2 + 1.5*b3*rss
    ga = 0.25-((pi*a*rs**5*b0)/24)*(rnd*rd-rn*rdd)/(rd*rd)
    g0 = (1./8.)*(z/special.i1(z))**2
    B = (9/16.)*ga - (3/64.)*(1-g0)-(16/15.)*A
    C = -0.75*ga + (9/16.)*(1-g0) - (16/5.)*A
    
    x2 = x*x
    x4 = x2*x2
    r = A*x4 + B*x2 + C
    if (x < 1e-6):
        r = r + (A*x4+(B+8/3.*A)*x2-C)
    elif (abs(2-x)>1e-6):
        r = r + (A*x4+(B+8/3.*A)*x2-C)*((4-x2)/(4*x))*log(abs((2+x)/(2-x)))
    return r

def gpv(kf, te, ge, igf=None, maxiter=100, gktol=0.05, sta=0.6,
        nz=200, zmin=1e-2, zmax=5e2, nw=2000, ipr=False):
    if igf is None:
        igf = interpgfun(-3, 3., 500, te, ge)
    de = kf**3/(3*pi**2)
    wp = sqrt(4*pi*de)
    ta = te*0.5*kf**2
    zs = 10**linspace(log10(zmin), log10(zmax), nz)
    wa = 5*(0.5*(ge+(ge**2+(pi*te)**2)**0.5))**0.5
    ks = zs*2*kf
    ks2 = ks*ks
    ks3 = ks2*ks
    aks = log(ks)
    ra = zeros((2,nz,nw))
    sk = zeros(nz)
    gk = zeros(nz)
    gk0 = gk.copy()
    za = zeros(nz)
    ws = zeros((nz,nw))
    for i in range(nz):
        wm = max(max(0.5*ks2[i],ks[i]*kf), wp)        
        wmin = wm-wa*ks[i]*kf
        wmax = wm+max(wa*ks[i]*kf,wp*3.0)
        if (wmin < 1e-6*wm):
            wmin = wm*1e-6
        ws[i] = linspace(wmin, wmax, nw)
        for j in range(nw):            
            ra[:,i,j] = rpa(zs[i], ws[i,j]/(ks[i]*kf), 0.0, kf, te, ge, igf)
    rc = ra.copy()

    niter = 0
    while (niter < maxiter):
        niter += 1
        im = elf((rc[0],rc[1]))
        for i in range(nz):
            k = ks[i]
            vk = 4*pi/(k*k)
            denom = de*pi*vk
            w = ws[i]
            dw = (w[1]-w[0])
            iy0 = im[i]/denom
            iy1 = exp(-w/ta)
            iy2 = 1.0/(1-iy1)
            iyp = iy0*iy2
            iy2 = iy1/(1-iy1)
            iym = iy0*iy2
            a = 1+rc[0,i]
            i0 = lastzero(a)
            ep0 = 0.0
            ep1 = 0.0
            em0 = 0.0
            em1 = 0.0
            w0 = 0.0
            if i0 < nw-1:
                i1 = i0+1
                j0 = max(0,i0-1)
                j1 = min(nw-1,i1+1)
                ak = median(array([(a[j1]-a[j0])/(w[j1]-w[j0]),
                                   (a[i1]-a[i0])/(w[i1]-w[i0]),
                                   (a[i1]-a[j0])/(w[i1]-w[j0]),
                                   (a[j1]-a[i0])/(w[j1]-w[i0])]))
                f = -a[i0]/(a[i1]-a[i0])
                w0 = w[i0]*(1-f) + w[i1]*f
                b0 = rc[1,i,i0]*(1-f) + rc[1,i,i1]*f        
                ak = abs(ak)
                ab = abs(b0)
                dw1 = ab/ak
                dw2 = dw1/dw
                if dw2 < 0.5:
                    dj = int(dw2+2)
                    j0 = max(0, i0-dj)
                    j1 = min(i1+dj, nw-1)
                    rkp = (iyp[j1]-iyp[j0])/(w[j1]-w[j0])
                    rkm = (iym[j1]-iym[j0])/(w[j1]-w[j0])
                    iyp[j0:j1] = iyp[j0] + (w[j0:j1]-w[j0])*rkp
                    iym[j0:j1] = iym[j0] + (w[j0:j1]-w[j0])*rkm
                    ep0 = integrate.simps(iyp[j0:j1], dx=dw)
                    em0 = integrate.simps(iym[j0:j1], dx=dw)
                    dws1 = w[j1]-w0
                    dws0 = w0-w[j0]
                    if (dw1 < 1e10*dws1):
                        dws1 = 1e10
                    else:
                        dws1 /= dw1
                    if (dw1 < 1e10*dws0):
                        dws0 = 1e10
                    else:
                        dws0 /= dw1
                    di = (arctan(dws1)+arctan(dws0))/ak
                    ewt = exp(-w0/ta)
                    ep1 = di/(1-ewt)/denom
                    em1 = di*ewt/(1-ewt)/denom
            sk[i] = integrate.simps(iyp, dx=dw)-ep0+ep1
            sk[i] += integrate.simps(iym, dx=dw)-em0+em1
        sk /= sk[-1]        
        for i in range(nz):
            iy = -(ks3*(sk-1.))/(4*pi*pi*de)
            ya = (ks2-ks2[i])**2/(4*ks*ks3[i])
            yb = log(abs(ks+ks[i])/(1e-31+abs(ks-ks[i])))            
            iy *= (5/6.-ks2/(2*ks2[i])+ya*yb)
            gk[i] = integrate.simps(iy, dx=aks[1]-aks[0])
        for j in range(nw):
            b = lfc((ra[0,:,j],ra[1,:,j]), (gk,za))
            rc[0,:,j] = b[0]
            rc[1,:,j] = b[1]
        dgk = sqrt(sum((gk-gk0)**2)/nw)
        if ipr:
            print('gpv iter: %d %g'%(niter, dgk))
        if (dgk < gktol):
            break
        if niter > 2:            
            gk = sta*gk + (1-sta)*gk0
        gk0[:] = gk
    return ks, ws, sk, gk, ra, rc

def vfn(r, zv, b2, ds, km):
    if (km <= 0):
        return -2*zv/r * exp(-r/ds)
    
    kmr = km*r
    if type(r) == float or type(r) == float64:
        if kmr < 1e-4:
            p = -2*zv*km*exp(-r/ds)
        else:
            p = -2*zv/r * exp(-r/ds)*(1-exp(-kmr))
    else:
        w = where(kmr < 1e-4)[0]
        p = zeros(len(r))
        if len(w) > 0:
            p[w] = -2*zv*km*exp(-r[w]/ds)
        w = where(kmr >= 1e-4)[0]
        if len(w) > 0:
            p[w] = -2*zv/r * exp(-r/ds)*(1-exp(-kmr[w]))
    return p

def cspfn(r, zv, b2, ds, km):
    r2 = r*r
    return vfn(r, zv, b2, ds, km) + b2/r2

def r2dvdr(r, zv, b2, ds, km):
    kd = 1/ds
    ex1 = exp(-kd*r)
    y = 1+kd*r
    if (km > 0):
        ex2 = exp(-km*r)
        y -= ex2*(1+(kd+km)*r)
    y *= zv*ex1
    y -= b2/r
    y *= 2
    return y

def r2dvdr2(r, zv, ds, km):
    kdr = r/ds
    p = -kdr**2 + kdr + 1.0
    if km > 0:
        kmr = km*r
        kmdr = kmr+kdr        
        p += exp(-kmr)*(kmdr**2-kmdr-1.0)
    return p

def csprp(zv, ds, km):
    r = 1.618033988749895*ds
    if km > 0:
        r0 = r
        r1 = r
        while(True):
            if (r2dvdr2(r1, zv, ds, km) < 0):
                break
            r1 *= 2.0
            
        while (r1-r0 > 1e-8*r):
            r = 0.5*(r0+r1)
            d = r2dvdr2(r, zv, ds, km)
            if d > 0:
                r0 = r
            elif d < 0:
                r1 = r
            else:
                break        
    return r

def csprz(zv, b2, ds, km, r0=0.0, r1=0.0, sr=-1):
    if zv < 0:
        return 0.0
    r = (sqrt(zv*zv+4*zv*b2/ds)-zv)/(2*zv/ds)
    if (r0 <= 0):
        if (r1 > 0):
            r0 = min(r, r1)
        else:
            r0 = r
        while(True):
            d = sr*r2dvdr(r0, zv, b2, ds, km)
            if (d < 0):
                break
            r0 /= 2.0
    if (r1 <= 0):
        r1 = max(r, r0)
        while(True):
            d = sr*r2dvdr(r1, zv, b2, ds, km)
            if (d > 0):
                break
            r1 *= 2.0

    while(r1-r0 > 1e-8*r):
        r = 0.5*(r0+r1)
        d = sr*r2dvdr(r, zv, b2, ds, km)
        if (d < 0):
            r0 = r
        elif (d > 0):
            r1 = r
        else:
            break
    return r

def csprm(zv, b2, ds, km, xt=1.0, r0=0.0, r1=0.0, sr=-1):
    r = (sqrt(zv*zv+b2*xt)-zv)/xt
    if r0 <= 0:
        if (r1 <= 0):
            r0 = r
        else:
            r0 = min(r,r1)
        while(True):
            d = sr*(xt-cspfn(r0, zv, b2, ds, km))
            if (d > 0):
                break
            r0 /= 2.0
    if r1 <= 0:
        r1 = max(r,r0)
        while(True):
            d = sr*(xt-cspfn(r1, zv, b2, ds, km))
            if (d < 0):
                break
            r1 *= 2.0
    while(r1-r0 > 1e-8*r):
        r = 0.5*(r0+r1)
        d = sr*(xt-cspfn(r, zv, b2, ds, km))
        if d > 0:
            r0 = r
        elif d < 0:
            r1 = r
        else:
            break

    return r

def cspbint(x, r0, b, rb2, b2, zv, ds, km):
    s = sin(x)
    c = cos(x)
    c2 = c*c
    rs = r0/cos(x)
    yd = (-c2+rb2*(1-vfn(rs,zv,b2,ds,km)))
    yi = s/sqrt(yd)
    return yi

def cspb(zv, b, ds, km, eps0=1e-2, eps1=1e-3):
    b2 = b*b
    if (zv > 0):
        rp = csprp(zv, ds, km)
        rz = csprz(zv, b2, ds, km, r0=rp, r1=0., sr=-1)
        d = cspfn(rz, zv, b2, ds, km)
        if (d >= 1.0):
            r0 = csprm(zv, b2, ds, km, xt=1.0, r0=rz, r1=0., sr=-1)
        else:
            r0 = csprm(zv, b2, ds, km, xt=1.0, r0=0., r1=0., sr=-1)
    else:
        rp = 0.0
        r0 = csprm(zv, b2, ds, km, xt=1.0, r0=0., r1=0., sr=-1)
    """
    rp = 0.0
    r0 = csprm(zv, b2, ds, km, xt=1.0, r0=0., r1=0., sr=-1)
    """
    rb = r0/b
    rb2 = rb*rb
    y = integrate.quad(cspbint, eps0, pi/2-eps1, args=(r0,b,rb2,b2,zv,ds,km))
    y = y[0]
    y += cspbint(eps0,r0,b,rb2,b2,zv,ds,km)*eps0
    y += cspbint(pi/2-eps1,r0,b,rb2,b2,zv,ds,km)*eps1
    #xs = linspace(1e-2, pi/2-1e-3, 10000)
    #ys = array([cspbint(x,r0,b,rb2,b2,zv,ds,km) for x in xs])
    #return 1-cos(pi-y*b*2),y,r0,xs,ys
    return 1-cos(pi-y*2)

def cspd(zv, b, ds, km):
    return cspb(zv,b,ds,0.0)-cspb(zv,b,ds,km)

def csp(v, zp, km, km1, dbx=1.25, nbmin=5, nbmax=50, tol=0.01):    
    zv = zp/v/v
    bs = array([5.0/km1])
    ds = 1e4*bs[0]
    ys = array([cspb(zv, bs[0], ds, 0.0)])
    ym = array([cspb(zv, bs[0], ds, km)])
    dx = log(dbx)
    n = 0
    tol0 = 0.1*tol;
    while (True):
        n = n+1
        b = bs[0]/dbx
        bs = append(b, bs)
        ys = append(cspb(zv, b, ds, km1), ys)
        ym = append(cspb(zv, b, ds, km), ym)
        if n > nbmax:
            break
        if (n > nbmin and
            ys[0] < tol0 and
            ym[0] < tol0):
            break
    n = 0
    while (True):
        n = n+1
        b = bs[-1]*dbx
        bs = append(bs, b)
        ys = append(ys, cspb(zv, b, ds, km1))
        ym = append(ym, cspb(zv, b, ds, km))
        if n > nbmax:
            break
        if n > nbmin and abs(ys[-1]-ym[-1]) < tol*ys[-1]:
            break
    bs2 = bs*bs
    yb = ys-ym
    y = integrate.simps(yb*bs2, dx=dx)*2*pi
    return y,bs,ys,yb

def dedxkm(v,kf,te,rs,zp=1.0,sc=1.0,se=1.0,si=0.0):
    azp = abs(zp)
    ste = sqrt(1+te/2)
    ve = kf*ste
    vr = sqrt(ve*ve + v*v)    
    vrp = vr/ve
    eta = (0.57*se*azp/vr)
    if si > 0:
        eta *= (vrp)**si
    #eta = 0.3*se*azp*rs/(vrp*ste)**si
    km1 = 2*vr
    km = km1/sqrt(1+eta*eta)
    km1 *= sc
    return array([km,km1])

def cspv(zp, de, te, vmin=0.1, vmax=1e2, nv=25):
    tu = te/consts.HARTREE_EV
    vx = exp(linspace(log(vmin), log(vmax), nv))
    vt = sqrt(tu)
    kf,d0,t0,ge,wp,rs,ga = spar(de, te)
    ef = 0.5*kf**2
    ds = sqrt(tu/wp**2)
    km = array([dedxkm(v, kf, t0, rs, zp=zp) for v in vx*vt])
    y = array([csp(vx[i]*vt, zp, km[i,0], km[i,1])[0] for i in range(nv)])
    y /= zp**2
    w = where(y < 1e-99)[0]
    if len(w) > 0:
        y[w] = 1e-99    
    return vx, y, ga, vt, ge, t0, ef, d0, ds

def cspintmu(u, v, fi, mu, vmin):
    v2 = v*v
    uv = u*v*mu
    vr = sqrt(u*u+v2-2*uv)
    w = where(vr < vmin)[0]
    if (len(w) > 0):
        vr[w] = vmin
    sig = exp(fi(log(vr)))
    im = vr*(u*mu - v)*sig
    return integrate.simps(im, dx=mu[1]-mu[0])
    
def dedxbc(v,vi,y,ga,vt,te,ge,ef,ds,
           nu=101,umax=25.0,umin=0.01,nmu=21):
    try:
        nv = len(v)
        sc = 0
    except:
        nv = 1
        v = array([v])
        sc = 1 
    r = zeros(nv)
    fi = interpolate.interp1d(log(vi), log(y), bounds_error=False,
                              fill_value='extrapolate')    
    u = linspace((umin*vt), (umax*vt), nu)
    vmin = 1e-3*umin*vt
    mu = linspace(-1., 1., nmu)
    xu = (0.5*u*u/ef-ge)/te
    xf = array([exp1p(t) for t in xu])*u*u
    for i in range(nv):        
        im = (array([cspintmu(t, v[i], fi, mu, vmin) for t in u]))
        iy = im*xf
        r[i] = -integrate.simps(iy, dx=u[1]-u[0])
    r /= 34.1893*ga**1.5*(te*ef)/ds
    if sc == 1:
        r = r[0]
    r *= (v/vt)**2
    return r

def barkas(v, kf, te, wp, zp=1.0, md=3):
    v2 = v*v
    y = 2*v2/wp
    if y < 1e6:
        y = 1/(1-exp(-1/y))
    y2 = y*y
    y8 = (1+8/y2)**0.5
    yx = log(0.25*(1+y8))
    yb = (5/3)*log(y) + (1/3)*(1+0.5*y2*(1-y8)+4*yx)
    
    v3 = v2*v
    v4 = v2*v2
    vp = wp/kf
    yb *= zp*wp*pi/v3
    if md > 0:
        ve = md*vp*sqrt(1+te)
        vf = v4/(ve**4+v4)
        yb *= vf
    
    return yb

def bloch(v, kf, te, wp, zp=1.0, md=3):
    y = zp/v
    y2 = y*y
    y4 = y2*y2
    yb = -y2*(1.202 - y2*(1.042-0.855*y2+0.343*y4))
    yb = min(0.0, yb)
    vp = wp/kf    
    if md > 0:
        ve = md*vp*sqrt(1+te)
        v4 = v**4
        yb *= v4/(ve**4+v4)
    return yb

def dedx1k(k,kmc,v,kf,wp,te,ge,igf,nw,igk=None,wi=0.0,zp=1.0):
    kkf = k*kf
    wm = max(max(0.5*k*k,kkf),wp)
    a0 = 5*((0.5*(ge+(ge**2+(pi*te)**2)**0.5))**0.5)*kkf
    w0 = wm - a0
    if (w0 < 0):
        w0 = 0.0
    w1 = wm + max(a0,3*wp)
    kv = k*v
    w1 = min(w1, kv)
    w0 = min(w0, 1e-1*w1)
    ws = linspace(w0,w1,nw)
    dw = ws[1]-ws[0]      
    u = ws/kkf
    z = k/(2.0*kf)
    r = array([rpa(z, t, wi/kkf, kf, te, ge, igf) for t in u])
    if (not igk is None):
        gk0 = giu(k, kf)
        if (type(igk) != type(1)):
            gk1 = igk(log(k))
            gka = cdiv((ws*gk1,wp*gk0),(ws, wp))
        else:
            gka = (gk0, 0.0)
        ra,rb = lfc((r[:,0],r[:,1]),gka)
        r[:,0] = ra
        r[:,1] = rb
    kma = [1.0, kmc]
    p = []
    for kmc in kma:
        r *= kmc
        ri = ws*elf((r[:,0],r[:,1]))
        a = 1+r[:,0]
        er0 = 0.0
        er1 = 0.0
        w0 = 0.0
        ab = 0.0
        ak = 0.0
        dws0 = 0.0
        dws1 = 0.0
        i0 = lastzero(a)
        if i0 < nw-1:
            i1 = i0+1
            j0 = max(0,i0-1)
            j1 = min(nw-1,i1+1)
            ak = median(array([(a[j1]-a[j0])/(ws[j1]-ws[j0]),
                               (a[i1]-a[i0])/(ws[i1]-ws[i0]),
                               (a[i1]-a[j0])/(ws[i1]-ws[j0]),
                               (a[j1]-a[i0])/(ws[j1]-ws[i0])]))
            f = -a[i0]/(a[i1]-a[i0])
            w0 = ws[i0]*(1-f) + ws[i1]*f
            b0 = r[i0,1]*(1-f) + r[i1,1]*f        
            ak = abs(ak)
            ab = abs(b0)
            dw1 = ab/ak
            dw2 = dw1/dw
            if dw2 < 0.5:
                dj = int(dw2+2)
                j0 = max(0, i0-dj)
                j1 = min(i1+dj, nw-1)
                rk = (ri[j1]-ri[j0])/(ws[j1]-ws[j0])
                ri[j0:j1] = ri[j0] + (ws[j0:j1]-ws[j0])*rk
                er0 = integrate.simps(ri[j0:j1], dx=dw)
                dws1 = ws[j1]-w0
                dws0 = w0-ws[j0]
                if (dw1 < 1e10*dws1):
                    dws1 = 1e10
                else:
                    dws1 /= dw1
                if (dw1 < 1e10*dws0):
                    dws0 = 1e10
                else:
                    dws0 /= dw1
                er1 = w0/ak*(arctan(dws1)+arctan(dws0))
        p.append(integrate.simps(ri, dx=dw) + er1-er0)
    #return p,ws,ri,r,er0,er1,ak,ab,w0,dws0,dws1
    return p

def linfun(x, a, b):
    return a*x+b

def fitkm(ks, ri):
    mr = max(ri)
    w = where(ri > 0.8*mr)[0]
    mr = median(ri[w])
    w = where(ri > 0.9*mr)[0]
    w0 = w[-1]
    w = where(ri > 0.1*mr)[0]
    w1 = w[-1]
    w = range(w0,w1+1)
    if (w1 <= w0+1):
        km = ks[w0]
        a = 10.0
    else:
        y = ri[w]/max(mr,max(ri[w])*1.01)
        y = log(y/(1-y))
        x = log(ks[w])
        p = optimize.curve_fit(linfun, x, y)[0]
        a = -p[0]
        km = exp(p[1]/a)
        a = max(2,a)
    return km,a

def dedx(ep,kf,de,te,ge,igf,kmin=0.25,kmax=2.5,nk=200,
         nw=400,igk=None,wi=0.0,zp=1.0):
    v2 = 2*(ep/(consts.HARTREE_EV*consts.AMU))
    v = sqrt(v2)
    wp = sqrt(4*pi*de)
    a0 = (0.5*(ge+(ge**2+(pi*te)**2)**0.5))**0.5
    k0 = wp/max(kf,v)
    k1 = 2*(v+a0*kf)
    ks = exp(linspace(log(kmin*k0), log(kmax*k1), nk))
    km = dedxkm(v, kf, te, (3/wp**2)**(1/3.), zp=zp)
    km2 = km**2
    ks2 = ks**2    
    kmc = (km2[0]/(km2[0]+ks2))/(km2[1]/(km2[1]+ks2))
    ri = array([dedx1k(ks[i],kmc[i],v,kf,wp,te,ge,igf,nw,
                       igk=igk,wi=wi,zp=zp) for i in range(nk)])
    dk = log(ks[1]/ks[0])
    r = integrate.simps(ri[:,0], dx=dk)
    rc = integrate.simps(ri[:,1], dx=dk)
    a = (2/pi)/wp**2
    r *= a
    rc *= a
    rb = max(0.0, barkas(v, kf, te, wp, zp=zp))
    # stopping number, with strong collision cutoff, barkas term.
    return r,rc,rb
    #return r,rc,ks,ri,km,kmc

# x = (kf, de, te, ge, z, igf, igk, wi, ep)
def dedx_loop(xs):
    if xs is None:
        return []
    r = []
    for x in xs:
        r.append(dedx(x[-1],x[0],x[1],x[2],x[3],x[5],
                      zp=x[4],igk=x[6],wi=x[7]))
    return r

# x = (kf, de, te, ge)
def igf_loop(xs):
    if xs is None:
        return []
    return [interpgfun(-3., 3., 500, x[2], x[3]) for x in xs]

# x = (kf, de, ge, te, igf)
def igk_loop(xs):
    if xs is None:
        return []
    igk = []
    for x in xs:
        ks,ws,sk,gk,ra,rc = gpv(x[0], x[2], x[3], maxiter=6, igf=x[4])
        r = interpolate.interp1d(log(ks), gk, bounds_error=False,
                                 fill_value=(gk[0],gk[-1]))
        igk.append(r)
    return igk

def dedxbc_loop(xs):
    if xs is None:
        return []
    rx = []
    for x in xs:
        vx,vy,ga,vt,ge,te,ef,d0,ds = cspv(x[0], x[1], x[2])
        r = dedxbc(x[-1], vx*vt, vy, ga, vt, te, ge, ef, ds)
        rx.append(r)
    return rx

def dist_work(a, np):
    n = len(a)
    b = [None]*np
    k = [None]*np
    for i in range(0, n, np):
        ni = min(i+np,n)
        for j in range(i,ni):
            if (k[j-i] is None):
                k[j-i] = []
                b[j-i] = []
            b[j-i].append(a[j])
            k[j-i].append(j)
    return k,b

def linlow(v, v2, y):
    z = y/v2
    k = argmax(z)
    if (k < 2):
        return y
    a = z[k-1]/v[k-1]
    z[:k-1] = v[:k-1]*a
    return z*v2

def dedxmp(es, ds, ts, mgk=0, wi=0.0, zp=1.0, np=16):
    try:
        ne = len(es)
    except:
        ne = 1
        es = [es]
    try:
        nd = len(ds)
    except:
        nd = 1
        ds = [ds]
    try:
        nt = len(ts)
    except:
        nt = 1
        ts = [ts]
    try:
        nz = len(zp)
    except:
        nz = 1
        zp = [zp]
    es = array(es)
    ds = array(ds)
    ts = array(ts)
    zp = array(zp)
    
    ndt = nd*nt*nz
    np1 = min(np, ndt)
    np2 = min(np, ndt*ne)
    p1 = Pool(processes=np1)
    p2 = Pool(processes=np2)
    t0 = time()
    tb = t0
    print('calculating gfun ...')
    xs = []
    ws = []
    vs = sqrt(es*2/(consts.HARTREE_EV*consts.AMU))
    for m in range(nz):
        for i in range(nd):
            for j in range(nt):
                kf,de,te,ge,wp,rs,ga = spar(ds[i], ts[j])
                xs.append((kf, de, te, ge, zp[m]))
                ws.append((zp[m],ds[i],ts[j],vs))
    k0, b0 = dist_work(xs, np)
    r0 = p1.map(igf_loop, b0)
    igfs = [None]*ndt
    for i in range(np):
        if k0[i] is None:
            continue
        for k in range(len(k0[i])):
            igfs[k0[i][k]] = r0[i][k]
    t1 = time()
    print('done %10.3E'%(t1-t0))
    t0 = t1
    igks = [None]*ndt
    print('calculating gpv ...')
    ys = []
    for i in range(ndt):
        ys.append(xs[i]+(igfs[i],))
    if mgk > 2:
        k1, b1 = dist_work(ys, np)
        r1 = p1.map(igk_loop, b1)
        for i in range(np):
            if k1[i] is None:
                continue
            for k in range(len(k1[i])):
                igks[k1[i][k]] = r1[i][k]
    t1 = time()
    print('done %10.3E'%(t1-t0))

    nc = ndt*ne
    zs0 = []
    zs1 = []
    zs2 = []
    for i in range(ndt):
        for j in range(ne):
            zs0.append(ys[i]+(None,wi,es[j]))
            zs1.append(ys[i]+(1,wi,es[j]))
            zs2.append(ys[i]+(igks[i],wi,es[j]))

    zs = (zs0,zs1,zs2)
    
    if mgk == 2:
        mgks = [0,1]
    elif mgk == 4:
        mgk = [0,1,2]
    else:
        mgks = [mgk]
    nms = len(mgks)
    rs = []
    rc = []
    rb = []
    t0 = t1
    print('calculating dedx ...')
    for imd in range(len(mgks)):
        md = mgks[imd]
        print('md=%d'%md)
        rs.append(zeros(nc))
        rc.append(zeros(nc))
        rb.append(zeros(nc))
        k2, b2 = dist_work(zs[md], np)
        r2 = p2.map(dedx_loop, b2)
        for i in range(np):
            if k2[i] is None:
                continue
            for k in range(len(k2[i])):
                r2a = r2[i][k]
                rs[imd][k2[i][k]] = r2a[0]
                rc[imd][k2[i][k]] = r2a[1]
                rb[imd][k2[i][k]] = r2a[2]
    t1 = time()
    print('done %10.3E'%(t1-t0))
              
    if len(mgks) > 1:
        rs = array(rs).reshape((nms,nz,nd,nt,ne))
        rc = array(rc).reshape((nms,nz,nd,nt,ne))
        rb = array(rb).reshape((nms,nz,nd,nt,ne))
    else:
        rs = rs[0].reshape((nz,nd,nt,ne))
        rc = rc[0].reshape((nz,nd,nt,ne))
        rb = rb[0].reshape((nz,nd,nt,ne))

    k3, b3 = dist_work(ws, np)
    r3 = p1.map(dedxbc_loop, b3)
    rx = [None]*ndt
    for i in range(np):
        if k3[i] is None:
            continue
        for k in range(len(k3[i])):
            rx[k3[i][k]] = r3[i][k]
    rx = array(rx).reshape((nz,nd,nt,ne))
    rs = rs.squeeze()
    rc = rc.squeeze()
    rb = rb.squeeze()
    rx = rx.squeeze()
    t1 = time()
    print('total time %10.3E'%(t1-tb))
    #rx is the strong binary collision correction term
    return rs,rc,rb,rx

"""
tabulate dedx. Columns:
Log(Te)
Log(Ne)
E(MeV)
L1, stopping number with no corrections.
L2 with LFC
L3 without LFC, with strong binary collision correction
L4 with LFC, and strong binary collision correction
L5 Barkas correction term.
L6 Strong binary collision correction term
"""
def savdedx(fn, ts, ds, es, rs, rc, rb, rx):
    v = sqrt(es)
    with open(fn, 'w') as f:
        nt = len(ts)
        nd = len(ds)
        ne = len(es)
        f.write('# nt = %2d\n'%nt)
        f.write('# nd = %2d\n'%nd)
        f.write('# ne = %2d\n'%ne)
        rs = rs.reshape((2,nd,nt,ne))
        rc = rc.reshape((2,nd,nt,ne))
        rb = rb.reshape((2,nd,nt,ne))
        rx = rx.reshape((nd,nt,ne))
        for i in range(nt):
            for j in range(nd):
                f.write('# T = %2d %12.5E\n'%(i, ts[i]))
                f.write('# D = %2d %12.5E\n'%(j, ds[j]))
                r0 = rs[0,j,i]
                r1 = rs[1,j,i]
                rc0 = rc[0,j,i]
                rc1 = rc[1,j,i]
                im = argmax(r0/es)
                rb0 = rb[0,j,i]
                rb1 = rb[1,j,i]
                rxi = rx[j,i]            
                ra0 = rc0+rxi
                ra1 = rc1+rxi
                for k in range(ne):
                    f.write('%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n'%(log10(ts[i]),log10(ds[j]),es[k],r0[k],r1[k],ra0[k],ra1[k],rb0[k],rxi[k]))
                f.flush()

def convdedx(fn0, fn1, bmd=0):
    r = loadtxt(fn0, unpack=1)
    r = r.reshape((9,30,1,50))
    ts = [10**r[0,0,0,0]]
    ds = 10**r[1,:,0,0]
    es = r[2,0,0]
    rs = zeros((2,30,1,50))
    rc = zeros((2,30,1,50))
    rb = zeros((2,30,1,50))
    rx = zeros((30,1,50))
    rs[0] = r[3]
    rs[1] = r[4]
    rc[0] = r[5]
    rc[1] = r[6]
    rb[0] = r[7]
    rb[1] = r[7]
    rx = r[8]
    if bmd > 0:
        v = sqrt(es*2e6/consts.HARTREE_EV/consts.AMU)
        for i in range(len(ds)):
            kf,de,te,ge,wp,ra,ga = spar(ds[i], ts[0])
            rb[0,i,0] = array([barkas(x, kf, te, wp, md=bmd) for x in v])
            rb[1,i,0] = rb[0,i,0]
    savdedx(fn1, ts, ds, es, rs, rc, rb, rx)
    
def tabdedx(fn, ts=[0.05], ds=[1e23], es=[1e8],
            tr=None, dr=None, er=None, wi=0.0, zp=1., np=16):
    if not tr is None:
        ts = 10**linspace(log10(tr[0]),log10(tr[1]),tr[2])
    else:
        ts = array(ts)
    if not dr is None:
        ds = 10**linspace(log10(dr[0]),log10(dr[1]),dr[2])
    else:
        ds = array(ds)
    if not er is None:
        es = 10**linspace(log10(er[0]),log10(er[1]),er[2])
    else:
        es = array(es)

    rs,rc,rb,rx = dedxmp(es*1e6, ds, ts, mgk=2, wi=wi, zp=zp, np=np)
    savdedx(fn, ts, ds, es, rs, rc, rb, rx)
    
def totdedx(r, d, t, m):
    xd = t[1,:,0]
    yd = t[m+3]
    ep = t[2,0]
    n = len(ep)
    v2 = 2*ep*1e6/(consts.HARTREE_EV*consts.AMU)
    cc = (4*pi/v2)*(consts.HARTREE_EV*(consts.RBOHR*1e-8)**2*1e15)
    di = log10(d/(4*pi*r**2*(consts.RBOHR*1e-8)**3))
    rx = zeros(n)
    dr = append(r[0],diff(r))
    nr = len(r)
    drx = zeros((n,nr))
    dri = zeros((n,nr))
    for i in range(n):
        w = where(yd[:,i] > 0)
        fi = interpolate.interp1d(xd[w], log(yd[w,i]),
                                  bounds_error=False,
                                  fill_value='extrapolate')
        dri[i] = exp(fi(di))
        drx[i] = dri[i]*dr*d*cc[i]

    rx = sum(drx, 1)
    return ep,rx,drx,dri

if __name__ == '__main__':
    te = float(sys.argv[2])
    wi = float(sys.argv[3])
    tabdedx(sys.argv[1], ts=[te],
            dr=[1e19, 1e29, 30],
            er = [1e-3, 1e2, 50],
            wi = wi,
            np = 20)
    
