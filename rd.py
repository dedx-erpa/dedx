from pylab import *
import matplotlib.pyplot as plt
from pfac import fac
from pfac import aa
from scipy import interpolate, integrate
import os, string
import dts

def rden(od):
    f = '%s/rho.functions'%od
    d = loadtxt(f, unpack=1, max_rows=1)
    zt = int(d[0])
    d = loadtxt(f, unpack=1, max_rows=1, skiprows=1)
    te = d[0]
    rho = d[1]
    d = loadtxt(f, unpack=2, max_rows=1, skiprows=2)
    n = int(d[0])
    rs = d[2]
    
    r = loadtxt(f, unpack=1, skiprows=3)
    s = r.shape
    n1 = int(s[0]*s[1]/n)
    r = transpose(r).reshape((n1,n))
    w = list(range(n1))
    w = w[:1]+w[2:]
    return r[w]

def raurho(ext='s'):
    d0 = loadtxt('output/z79s/AuRho1.txt', unpack=1)
    d1 = loadtxt('output/z79s/AuRho2.txt', unpack=1)
    d = rden('output/z79%s'%ext)
    r = d1[0][0]
    rs = d1[0][-1]
    rs0 = d[0][-1]
    a = rs0/rs    
    f0 = interpolate.interp1d(d0[0]*a, d0[1], bounds_error=False,
                              fill_value='extrapolate')
    f1 = interpolate.interp1d(d1[0]*a, d1[1], bounds_error=False,
                              fill_value='extrapolate')
    w0 = where((d[0] >= 0.015) & (d[0] < r))
    w1 = where(d[0] >= r)
    y = d[1].copy()
    y[w0] = f0(d[0][w0])
    y[w1] = f1(d[0][w1])
    fy = interpolate.interp1d(d[0], y, bounds_error=False, fill_value='extrapolate')
    b = integrate.quad(fy, 0.0, rs0)
    b = 79./b[0]
    y *= b
    return array((d[0],d[1],d[2],y)),d0,d1,d

def intrho(r, d):
    fi = interpolate.interp1d(r, d, bounds_error=False,
                              fill_value='extrapolate')
    return integrate.quad(fi, r[0], r[-1])

def waurho(idir, odir):
    r,d0,d1,d = raurho(idir)
    z = 79
    te = 0.05
    de = 12.56
    rs = r[0][-1]
    nr = len(r[0])
    rmin = r[0][0]
    h = np.log(rs/rmin)/(nr-1)
    xr = np.arange(nr)
    rr = r[0]
    da = r[3]
    df = r[2]

    ofn = 'output/z79%s/rho.functions'%odir
    with open(ofn, 'w') as f:
        f.write('%5d %6d %6d\n'%(z, 1, 1))
        f.write('%12.5E %12.5E\n'%(te, de))
        f.write('%4d %12.5E %12.5E\n'%(nr, h, rs))
        for i in range(nr):
            f.write('%13.7E '%rr[i])
            if (i+1)%5 == 0:
                f.write('\n')
        for i in range(nr):
            f.write('%13.7E '%1.0)
            if (i+1)%5 == 0:
                f.write('\n')
        for i in range(nr):
            f.write('%13.7E '%da[i])
            if (i+1)%5 == 0:
                f.write('\n')
        for i in range(nr):
            f.write('%13.7E '%df[i])
            if (i+1)%5 == 0:
                f.write('\n')
                    
def rdedx(od, header=None):
    f = '%s/dedx.dat'%od
    nzt = int(loadtxt(f, unpack=1, comments='A', max_rows=1, usecols=3))
    zt = loadtxt(f, comments='A', max_rows=1, skiprows=1,
                usecols=range(3,3+nzt), ndmin=1).flatten()
    wt = loadtxt(f, comments='A', max_rows=1, skiprows=2,
                usecols=range(3,3+nzt), ndmin=1).flatten()
    d = loadtxt(f, comments='A', max_rows=7, skiprows=3, usecols=3, unpack=1)
    s = loadtxt(f, unpack=1, comments='A', max_rows=7, skiprows=3,
                usecols=1, dtype=str)
    h = {'nzt':nzt, 'zt':zt, 'wt':wt}
    for i in range(len(d)):
        if s[i] in ['mep', 'mloss']:
            h[s[i]] = int(d[i])
        else:
            h[s[i]] = d[i]
    if not header is None:
        hs = header.strip()
        if hs == '':
            return h
        return h[hs]
    
    z = h['zp']
    r = loadtxt(f, unpack=1)
    #barkas effect
    i = argmax(r[1])
    ek = r[0]*1e3
    v2 = (2*ek*1e3/27.21/1836.)
    cc = 4*pi/v2 * z * 27.21*0.53e-8**2*1e15
    cek = (ek/(ek+ek[i]))**2
    #cl0 = 0.001*ek*(ek/(ek+ek[i]))**2
    #cl1 = 1.5/ek**0.4 + 4.5e4/(z*ek**1.6)
    cl0 = 0.00182*ek*cek
    cl1 = 19.6/ek
    cl01 = cc*cl0*cl1/(cl0+cl1) 
    #r[1] += cl01
    #r[2] += cl01
    return append(r,cl01.reshape((1,len(cl0))),axis=0)

def rcdedx(od):
    d = rden(od)
    f = '%s/rdedx.dat'%od    
    r = loadtxt(f, unpack=1, usecols=(2,3,4))
    nd = len(d[0])
    nx = int(len(r[0])/nd)
    r = r.reshape((3,nx,nd))
    return d,r

def plot_rcdedx(d,r,i,nm=0):
    clf()
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    if nm == 0:
        ym = 1.0
    else:
        ym = max(r[1,i]/d[0])
    plot(d[0], ym*d[1]/max(d[1]), label=r'$4\pi r^2 \rho(r)$')
    y = r[1,i]/d[0]/d[1]
    plot(d[0], ym*y/max(y), label=r'$L(r)$')
    y = r[1,i]/d[0]
    plot(d[0], ym*y/max(y), label=r'$4\pi r^2 \rho(r) L(r)$',marker='.')
    plot(d[0], ym*r[2,i], label=r'$\int_0^r 4\pi r^2\rho(r)L(r)dr$')
    xlabel('radius (atomic unit)')
    ylabel('arb. unit')
    title('E=%4.2f keV'%(1e3*r[0,i,0]))    
    legend()

def prange(zs=[6,13,47,79]):
    clf()
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    iz = 0
    if type(zs) != type([]):
        zs = [zs]
        iz = 1
    for z in zs:
        if type(z) == type(0):
            a = fac.ATOMICSYMBOL[z]
        else:
            a = z
            z = 0
        od = 'data/%s'%a
        r = rdedx(od)
        loglog(r[0], r[2], label=a)
        if os.path.exists('%s/pstar.txt'%od):
            d = loadtxt('%s/pstar.txt'%od,
                        unpack=1, skiprows=8)
            w = where((d[0] > r[0][0])&(d[0] < r[0][-1]))[0]
            loglog(d[0][w], d[4][w]*1e3, marker='o',
                   fillstyle='none',           
                   markeredgewidth=0.75, markersize=2.5,
                   linestyle='none', label=a+' PSTAR')
    legend()
    xlabel('Energy (MeV)')
    ylabel('Proton Range (mg/cm$^2$)')
    if iz > 0:
        savefig('%s/range.pdf'%od)
    else:
        savefig('range.pdf')
    
def rpstar(od):
    f = '%s/rho.functions'%od
    d = loadtxt(f, unpack=1, max_rows=1)
    zt = int(d[0])
    
    f = '%s/pstar.txt'%od

    r = loadtxt(f, unpack=1, skiprows=8)
    r[[1,2,3]] *= 1e6/(1./(fac.ATOMICMASS[zt]*1.67e-24))*1e15
    return r

def e2v(e):
    return sqrt(2*e*1e6/27.21/1836.0)

def tofloat(s):
    try:
        return float(str(s).replace(',',' ').strip())
    except:
        return 0.0

def int_range(r, dx=0.01, m=0.0, z=0):
    x = log10(r[0])
    y = r[1].copy()
    if z > 0 and m == 0.0:
        m = fac.ATOMICMASS[z]
    if m > 0:
        y = y*1e-21/(m*1.67e-21)
    x = append(x[0]-3.0, x)
    y = append(10**(-3.0*log10(y[1]/y[0])/(x[2]-x[1]))*y[0],y)
    x0 = x[0]
    x1 = x[-1]
    xs = arange(x0, x1+dx*0.1, dx)
    f = interpolate.interp1d(x, log(y))
    ys = 10**xs/exp(f(xs))
    n = len(xs)
    rs = zeros(n)
    h = dx*log(10)*0.5
    for i in range(1, n):
        rs[i] = rs[i-1]+h*(ys[i-1]+ys[i])

    f = interpolate.interp1d(xs[1:], log(rs[1:]), bounds_error=False, fill_value='extrapolate')
    return exp(f(x[1:]))

def eloss(x, y, e, dx):
    xx = log(x)
    yy = log(y)
    f0 = interpolate.interp1d(xx, yy, bounds_error=False, fill_value='extrapolate')
    f1 = interpolate.interp1d(yy, xx, bounds_error=False, fill_value='extrapolate')
    dx1 = exp(f0(log(e)))-dx
    if dx1 <= 0:
        return e
    return e - exp(f1(log(dx1)))

def plot_eloss(a, dx, ds=['dt0', 'dt1', 'dt2'],
               labs=['cold', 'low current', 'high current'], mp=2, e0=0.8, e1=1.2):
    x = linspace(e0, e1, 25)
    clf()
    try:
        z = fac.ATOMICSYMBOL.index(a)
        c = 1e-18/(fac.ATOMICMASS[z]*1.67e-21)
    except:
        if a == 'Mylar':
            c = 1e-18/(((4+5*12+2*16)/11)*1.67e-21)
        else:
            c = 0.0
    cols = ['k','b','r','m','g','y']
    for i in range(len(ds)):
        r = rdedx('tmp/%s%s'%(ds[i],a))
        y = array([eloss(r[0]*mp, r[2]*mp, e, dx) for e in x])
        plot(x, y*1e3, label=labs[i], marker='o', color=cols[i])
    if mp == 2:
        xlabel('Deuteron Energy (MeV)')
    else:
        xlabel('Proton Energy (MeV)')
    ylabel('Energy Loss (keV)')

    title(r'%s $\rho$dx=%g'%(a,dx))
    legend()
    savefig('eloss%d_%s.png'%(mp,a))
    savefig('eloss%d_%s.pdf'%(mp,a))
    
def cmp_dedx(zt, xsc='e', dpass=1, atima=1, pstar=1, rpa=1, orpa=1):
    clf()    
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    if type(zt) == type(0):
        a = fac.ATOMICSYMBOL[zt]
    else:
        a = zt
        try:
            zt = fac.ATOMICSYMBOL.index(a)
        except:
            zt = 0
    d = 'data/%s'%a
    lab = 'eRPA'

    d0 = d
    if (os.path.exists(d0+'/atima.txt') and atima==1):
        r = loadtxt(d0+'/atima.txt', unpack=1)
        x = r[0]
        if (xsc == 'v'):
            x = e2v(x)
            plot(x, r[1], label='ATIMA')
        else:
            semilogx(x, r[1], label='ATIMA')
    if (os.path.exists(d0+'/dpass.txt') and dpass==1):
        r = loadtxt(d0+'/dpass.txt', unpack=1, skiprows=13)
        x = r[0]
        if (xsc == 'v'):
            x = e2v(x)
            plot(x, r[1], label='DPASS')
        else:
            semilogx(x, r[1], label='DPASS')
        
    if (xsc == 'e'):
        xlabel('Energy (MeV)')
        xlim(1e-3,1e2)
    else:
        xlabel('Velocity (a.u.)')
        xlim(0, 20.0)

    if (os.path.exists(d0+'/pstar.txt') and pstar==1):
        r = rpstar(d0)
        x = r[0]
        if (xsc == 'v'):
            x = e2v(x)
            plot(x, r[1], marker='o', markerfacecolor='none', color='k', label='PSTAR')
        else:
            semilogx(x, r[1], marker='o', markerfacecolor='none', color='k', label='PSTAR')
            
    fds = []
    fds.append(a+'.txt')
    fds.append(a.lower()+'.txt')
    fds.append(a.upper()+'.txt')
    fds.append('h'+a.lower()+'.txt')
    fds.append('H'+a+'.txt')
    fds.append('H'+a.upper()+'.txt')
    if a == 'H4C5O2':
        fds.append('MYLAR.txt')
    for fd in fds:
        try:
            if zt == 40:
                r = loadtxt('data/iaea/'+fd,
                            unpack=1,usecols=(1,2),dtype=str)
            else:
                r = loadtxt('data/iaea/'+fd,
                            unpack=1,usecols=(0,1),dtype=str)
            x = array([tofloat(r[0][i]) for i in range(len(r[0]))])
            y = array([tofloat(r[1][i]) for i in range(len(r[1]))])
            i = where((x > 0)&(y > 0))
            if ((zt != 14 and x[0] == 0 and y[0] == 0 and x[1] == 0) or
                max(y) < 0.1):
                y = y*fac.ATOMICMASS[zt]*1.67
            else:
                x = x/1e3
            plot(x, y, marker='o', color='r',
                 markersize=3, markeredgewidth=0.5,
                 fillstyle='none', linestyle='none', label='Exp.')
            break
        except:
            pass
                    
    y = rdedx(d)
    x = y[0]
    if (xsc == 'v'):
        x = e2v(x)
        plot(x, y[1], linewidth=2, label=lab)
    else:
        semilogx(x, y[1], linewidth=2, label=lab)
    if (os.path.exists(d+'0/dedx.dat') and rpa==1):
        y = rdedx(d+'0')
        x = y[0]
        if (xsc == 'v'):
            x = e2v(x)
            plot(x, y[1], linestyle='--', linewidth=3, label='RPA')
        else:
            semilogx(x, y[1], linestyle='--', linewidth=3, label='RPA')
    #if (os.path.exists(d+'0/dedx.ppd')):
    #    y = loadtxt(d+'0/dedx.ppd', unpack=1)
    #    y[1] *= 1e15
    #    plot(y[0], y[1], linestyle='dashdot', linewidth=3, label='RPA+QEOS')
    if (os.path.exists(d+'/odedx.dat') and orpa==1):
        y = loadtxt(d+'/odedx.dat',unpack=1)
        x = y[0]
        if (max(x) > 1000):
            x /= 1e3
        plot(x, y[1], linestyle='dashdot', linewidth=3, label='Wang et al')

    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    title('Target=%s'%a)    
    legend()

def plot_miter(zs = dts.ss):
    a = aa.AA()
    yz = []
    for z in zs:
        if type(z) == type(0):
            z = fac.ATOMICSYMBOL[z]
        od = 'data/%s'%z
        h = rdedx(od, header='')
        zt = h['zt']
        for x in zt:
            r = a.rpot('%s/%s'%(od,fac.ATOMICSYMBOL[int(x)]),header='')
            yz.append(r['miter'])
    xz = arange(len(yz))
    clf()
    plot(xz, yz, marker='.')
    
def gen_plots(zs = dts.ss, mrc=0):
    for z in zs:
        print(z)
        prange(z)
        cmp_dedx(z, atima=0, dpass=0)
        if type(z) == type(0):
            od = 'data/%s'%(fac.ATOMICSYMBOL[z])
        else:
            od = 'data/%s'%z
        savefig('%s/dedx.png'%od)
        savefig('%s/dedx.pdf'%od)
        if mrc > 0:
            d,r = rcdedx(od)        
            for m in range(10):
                plot_rcdedx(d, r, m)
                savefig('%s/rdedx%d.png'%(od,m))
                savefig('%s/rdedx%d.pdf'%(od,m))
        pause(0.25)
        
def pden0(z):
    od = 'data/%s'%fac.ATOMICSYMBOL[z]
    d = rden(od)
    clf()
    plot(d[0], d[1])
    if (z == 79):
        r = loadtxt('../dedx.old/output/z79s/au_rho.txt', unpack=1)
        plot(r[0], r[1], marker='.')
    xlabel('radius (a.u.)')
    ylabel(r'$4\pi r^2\rho(r)$ (a.u.)')

def pden(s, rs=1, ys=1, fd=0, dd='data', op=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    od = '%s/%s'%(dd,s)
    h = rdedx(od, header='')
    zt = h['zt']

    if op == 0:
        clf()
    for iz in range(len(zt)):
        z = int(zt[iz])
        if (h['nzt'] == 1):
            r = rden(od)
        else:
            r = rden(od+'/'+fac.ATOMICSYMBOL[z])
        x = r[0].copy()
        d = r[1].copy()
        f = r[2].copy()
        if ys == 1:
            d = d/(4*pi*x**2)
            f = f/(4*pi*x**2)
        if rs == 1:
            x /= x[-1]
        elif rs == 2:
            x = sqrt(x/x[-1])
        plot(x, d, label=fac.ATOMICSYMBOL[z])
        if fd > 0:
            plot(x, f, label=fac.ATOMICSYMBOL[z]+' free')
        if (ys == 1):
            yscale('log')

    legend()
    if rs == 1:
        xlabel('r/rs')
    elif rs == 2:
        xlabel(r'$\sqrt{r/rs}$')
    else:
        xlabel('radius')
    if ys == 1:
        ylabel(r'$\rho(r)$')
    else:
        ylabel(r'$4\pi r^2 \rho(r)$')
    title(s)
    
def gen_den_plots():    
    zs = [3, 4, 5, 6, 7, 8, 10, 13, 14, 18, 22,
          26, 28, 29, 32, 36, 42, 47,
          50, 54, 64, 74, 78, 79, 82, 92]
    for z in zs:
        pden(z)
        savefig('data/%s/rho.pdf'%(fac.ATOMICSYMBOL[z]))
    
def pcomp(c):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    r0 = loadtxt('tc%s/dedx.dat'%c, unpack=1)
    r1 = loadtxt('th%s/dedx.dat'%c, unpack=1)
    if c == 'BN':
        m = (10.31*0.5+14.01*0.5)*1.67e-24
        d = None
    elif c == 'Al2O3':
        m = (27.*0.4 + 16*0.6)*1.67e-24
        d = loadtxt('pstar_al2o3.txt', unpack=1, skiprows=8)
    elif c == 'SiO2':
        m = (28.1 + 16.0*2)/3*1.67e-24
        d = loadtxt('pstar_sio2.txt', unpack=1, skiprows=8)
    clf()
    loglog(r0[0], r0[1]*1e-21/m, marker='.', label='Cold')
    loglog(r1[0], r1[1]*1e-21/m, marker='.', label='T=1 keV')
    if not d is None:
        loglog(d[0], d[1], label='PSTAR')
    legend()
    xlim(5e-4,2e2)
    ylim(1, 1e3)
    xlabel('Energy (MeV)')
    ylabel(r'dE/dx (MeV cm$^2$/g)')
    if c == 'BN':
        title('Boron Nitride')
        savefig('BN.png')
        savefig('BN.pdf')
    elif c == 'Al2O3':        
        title('Aluminum Oxide')
        savefig('Al2O3.png')
        savefig('Al2O3.pdf')
    elif c == 'SiO2':        
        title('Silicon Dioxide')
        savefig('SiO2.png')
        savefig('SiO2.pdf')

def pbd():
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    nd = 6
    ds = 4*10**linspace(0.0, 5.0, nd)
    clf()
    for i in range(nd):
        r = rdedx('tmp/d%dB'%i)
        loglog(r[0], r[1], label=r'$\rho$=%d g/cc'%(int(ds[i])))
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')

    legend()
    
    savefig('tmp/bd_dedx.pdf')
    
def pbt(ide):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.234, 2.34]
    ts = [0.025, 5.0, 10., 15., 20., 25., 30., 100., 300., 1e3, 3e3, 10e3]
    nd = len(ds)
    nt = len(ts)
    clf()
    for i in range(nt):
        if i > 0 and i < 6:
            continue
        r = rdedx('tmp/dtB/d%dt%02dB'%(ide,i))
        loglog(r[0], r[1], label=r'T=%g eV'%(int(ts[i])))
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    title(r'$\rho$=%g g/cc'%ds[ide])
    legend()
    
    savefig('tmp/bt%d_dedx.pdf'%(ide))
    
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    clf()
    for i in range(nt):
        r = rdedx('tmp/dtB/d%dt%02dB'%(ide,i))
        loglog(r[0], r[2]/1e3, label=r'T=%g eV'%(int(ts[i])))
    xlabel('Energy (MeV)')         
    ylabel('Proton Range (g/cm$^2$)')          
    title(r'$\rho$=%g g/cc'%ds[ide])
    legend()
    
    savefig('tmp/bt%d_range.pdf'%ide)

def pct(ide=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.05, 0.1, 0.5]
    ts = [0.025, 10, 20, 30]
    nd = len(ds)
    nt = len(ts)
    clf()
    a = 1e-18/(12*1.67e-24*1e6)
    for i in range(nt):
        r = rdedx('tmp/d%dt%02dC'%(ide,i))
        plot(r[0], r[1]*a, label=r'T=%g eV'%(int(ts[i])))
    r0 = rdedx('data/C')
    r1 = rpstar('data/C')
    plot(r0[0], r0[1]*a, label='Solid')
    plot(r1[0], r1[1]*a, label='PSTAR')
    xlim(0, 1.0)
    ylim(0, 1.5)
    d = loadtxt('malko_fig.txt', unpack=1)
    plot(d[0][::2]/1e3, d[1][::2], marker='o', color='k', label='DFT', linestyle='none')
    errorbar(d[0][::2]/1e3, d[1][::2], yerr=(d[1][1::2]-d[1][::2]), capsize=3, marker='o', color='k', linestyle='none')
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx (keV/($\mu$g/cm$^2$))')
    title(r'$\rho$=%g g/cc'%ds[ide])
    legend()
    
    savefig('tmp/ct%d_dedx.png'%(ide))
    savefig('tmp/ct%d_dedx.pdf'%(ide))
    
def paut(ide=0, nde=0, md=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.19311, 1.9311, 19.311]
    ts = 10**linspace(-1,4.,35)

    nd = len(ds)
    if nde == 0:
        nde = nd
    nt = len(ts)
    clf()
    rt = zeros((nd,nt))
    ry = zeros((nd,nt))
    rz = zeros((nd,nt))
    rf = zeros((nd,nt))
    mi = zeros((nd,nt))
    a = aa.AA()
    for j in range(nd):
        for i in range(nt):
            r = rdedx('tmp/dtAu/d%dt%02dAu'%(j,i))
            fi = interpolate.interp1d(r[0], r[2])
            ry[j,i] = fi(1.6)
            h = a.rden('tmp/dtAu/d%dt%02dAu/Au'%(j,i), header='')
            rt[j,i] = h['T']
            rz[j,i] = h['ub']
            rf[j,i] = h['zf']
            h = a.rpot('tmp/dtAu/d%dt%02dAu/Au'%(j,i), header='')
            mi[j,i] = h['miter']
            if (mi[j,i]>100):
                print([j,i,mi[j,i]])
    for j in range(ide,nde+ide):
        lab = r'$\rho=%g$'%ds[j]
        if md == 0:
            plot(rf[j], ry[j], label=lab, marker='o')
        elif md == 1:
            semilogx(rt[j], rz[j], label=lab, marker='o')
        elif md == 2:
            semilogx(rt[j], rf[j], label=lab, marker='o')
        else:
            semilogx(rt[j], mi[j], label=lab, marker='o')
    if md == 0:
        xlabel('Zbar')                   
        ylabel(r'range (mg/cm$^2$)')
        xlim(-1,75)
        ylim(10,50)
    elif md == 1:
        xlabel('Temerature')
        ylabel(r'$\mu$')              
    elif md == 2:
        xlabel('Temperature')
        ylabel('Zbar')
    else:
        xlabel('Temperature')
        ylabel('max iter')
    title('Proton in Au')
    legend()
    
    savefig('tmp/aut%d_range.pdf'%(ide))    
    
def pnit(ide=0, nde=0, md=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.0891, 0.891, 8.91]
    ts = 10**linspace(-1,4.,35)

    nd = len(ds)
    if nde == 0:
        nde = nd
    nt = len(ts)
    clf()
    ry = zeros((nd,nt))
    rz = zeros((nd,nt))
    rf = zeros((nd,nt))
    rt = zeros((nd,nt))
    mi = zeros((nd,nt))
    a = aa.AA()
    for j in range(nd):
        for i in range(nt):
            r = rdedx('tmp/dtNi/d%dt%02dNi'%(j,i))
            fi = interpolate.interp1d(r[0], r[2])
            ry[j,i] = fi(1.6)
            h = a.rden('tmp/dtNi/d%dt%02dNi/Ni'%(j,i), header='')
            rt[j,i] = h['T']
            rz[j,i] = h['ub']
            rf[j,i] = h['zf']
            h = a.rpot('tmp/dtNi/d%dt%02dNi/Ni'%(j,i), header='')
            mi[j,i] = h['miter']
            if (mi[j,i]>100):
                print([j,i,mi[j,i]])
    for j in range(ide,nde+ide):
        lab = r'$\rho=%g$'%ds[j]
        if md == 0:
            plot(rf[j], ry[j], label=lab, marker='o')
        elif md == 1:
            semilogx(rt[j], rz[j], label=lab, marker='o')
        elif md == 2:
            semilogx(rt[j], rf[j], label=lab, marker='o')
        else:
            semilogx(rt[j], mi[j], label=lab, marker='o')
    if md == 0:
        xlabel('Zbar')                   
        ylabel(r'range (mg/cm$^2$)')
        xlim(-1,26)
        ylim(3,20)
    elif md == 1:
        xlabel('Temerature')
        ylabel(r'$\mu$')              
    elif md == 2:
        xlabel('Temperature')
        ylabel('Zbar')
    else:
        xlabel('Temperature')
        ylabel('max iter')
        
    title('Proton in Ni')
    legend()
    
    savefig('tmp/nit%d_range.pdf'%(ide))
    
def palt(ide=0, nde=0, md=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.027,0.27, 2.7]
    ts = 10**linspace(-1,4.0,35)

    nd = len(ds)
    nt = len(ts)
    if nde == 0:
        nde = nd
    clf()
    rt = zeros((nd,nt))
    ry = zeros((nd,nt))
    rz = zeros((nd,nt))
    rf = zeros((nd,nt))
    mi = zeros((nd,nt))
    a = aa.AA()
    for j in range(nd):
        for i in range(nt):
            r = rdedx('tmp/dtAl/d%dt%02dAl'%(j,i))
            fi = interpolate.interp1d(r[0], r[2])
            ry[j,i] = fi(1.6)
            h = a.rden('tmp/dtAl/d%dt%02dAl/Al'%(j,i), header='')
            rt[j,i] = h['T']
            rz[j,i] = h['ub']
            rf[j,i] = h['zf']
            h = a.rpot('tmp/dtAl/d%dt%02dAl/Al'%(j,i), header='')
            mi[j,i] = h['miter']
            if (mi[j,i]>100):
                print([j,i,mi[j,i]])
    for j in range(ide,nde+ide):
        lab = r'$\rho=%g$'%ds[j]
        if md == 0:
            plot(rf[j], ry[j], label=lab, marker='o')
        elif md == 1:
            semilogx(rt[j], rz[j], label=lab, marker='o')
        elif md == 2:
            semilogx(rt[j], rf[j], label=lab, marker='o')
        else:
            semilogx(rt[j], mi[j], label=lab, marker='o')
    if md == 0:
        xlabel('Zbar')                   
        ylabel(r'range (mg/cm$^2$)')
        xlim(-1,12)
        ylim(3,12)
    elif md == 1:
        xlabel('Temerature')
        ylabel(r'$\mu$')              
    elif md == 2:
        xlabel('Temperature')
        ylabel('Zbar')
    else:
        xlabel('Temperature')
        ylabel('max iter')
        
    title('Proton in Al')
    legend()
    
    savefig('tmp/alt%d_range.pdf'%(ide))

def all_figs():
    d,r = rcdedx('data/Al')
    for i in [0,5,9]:
        plot_rcdedx(d, r, i)
        savefig('rc%02d.pdf'%i)
    pbd()
    pbt(0)
    pbt(1)
    pct()
    paut()
    pnit()
    palt()
    
def pfit():
    clf()
    zs = [6,13,14,29,47,64,79]
    cs = ['k','r','g','b','y','c','m']
    for i in range(len(zs)):
        d = 'data/%s'%fac.ATOMICSYMBOL[zs[i]]
        r = rdedx(d)
        p = rpstar(d)
        yp = interp(r[0], p[0], p[1])
        semilogx(r[0], yp, marker='o', markersize=3, color=cs[i],
                 markerfacecolor='none', linestyle='none')
        semilogx(r[0], r[1], color=cs[i], label=fac.ATOMICSYMBOL[zs[i]])
    legend()
    
def diffpr(f0, f1):
    a = aa.AA()
    r0 = a.rden(f0, header='')
    r1 = a.rden(f1, header='')
    de = (r1['vf']-r0['vf'])
    dv = (4*pi/3)*(r1['rps']**3-r0['rps']**3)
    c = 1.6e-19/(0.53e-10**3)/1e9
    r = c*(-de/dv)
    v = c*(((0.5*(r0['ve']+r0['vn'])+2*r0['ek'])/3-0.5*r0['vx'])/((4*pi/3)*r0['rps']**3))
    return de,dv,r,v
