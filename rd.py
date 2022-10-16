from pylab import *
import matplotlib.pyplot as plt
from pfac import fac
from pfac import aa
from scipy import interpolate, integrate
import os, string
import dts

def zwc(s):
    a = fac.ATOMICSYMBOL
    s = s+'X'
    n = len(s)
    zs = []
    ws = []
    i = 0
    j = 0
    for k in range(i+1,n):
        if j == i and s[k].isdigit():
            j = k
            continue
        if s[k].isupper():
            if j == i:
                zs.append(a.index(s[i:k]))
                ws.append(1.0)
            else:
                zs.append(a.index(s[i:j]))
                ws.append(float(s[j:k]))
            i = k
            j = k
    return zs,ws

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

    return r[[0,2,3,4]]

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

def plot_rcdedx(d,r,i):
    clf()
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    plot(d[0], d[1]/max(d[1]), label=r'$4\pi r^2 \rho(r)$')
    y = r[1,i]/d[0]/d[1]
    plot(d[0], y/max(y), label=r'$L(r)$')
    y = r[1,i]/d[0]
    plot(d[0], y/max(y), label=r'$4\pi r^2 \rho(r) L(r)$')
    plot(d[0], r[2,i], label=r'$\int_0^r 4\pi r^2\rho(r)L(r)dr$')
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
        r = rdedx('%s%s'%(ds[i],a))
        y = array([eloss(r[0]*mp, r[2]*mp, e, dx) for e in x])
        plot(x, y*1e3, label=labs[i], marker='o', color=cols[i])
        if c > 0:
            f = interpolate.interp1d(r[0]*mp,r[1]*c)
            y1 = f(x)*dx
            for j in range(10):
                y1 = f(x-0.5*y1/1e3)*dx
            plot(x, y1, linestyle='--', color=cols[i])
    if mp == 2:
        xlabel('Deuteron Energy (MeV)')
    else:
        xlabel('Proton Energy (MeV)')
    ylabel('Energy Loss (keV)')

    title(r'%s $\rho$dx=%g'%(a,dx))
    legend()
    savefig('eloss%d_%s.png'%(mp,a))
    savefig('eloss%d_%s.pdf'%(mp,a))
    
def cmp_dedx(zt, xsc='e', dpass=1, atima=1, pstar=1):
    clf()    
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    if type(zt) == type(0):
        a = fac.ATOMICSYMBOL[zt]
    else:
        a = zt
        zt = 0
    d = 'data/%s'%a
    lab = 'eRPA'
    y = rdedx(d)
    x = y[0]
    if (xsc == 'v'):
        x = e2v(x)
        plot(x, y[1], label=lab)
    else:
        semilogx(x, y[1], label=lab)

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
    for fd in fds:
        try:
            if zt == 40:
                r = loadtxt('/home/mfgu/dd/dedx_notes/iaea/'+fd,
                            unpack=1,usecols=(1,2),dtype=str)
            else:
                r = loadtxt('/home/mfgu/dd/dedx_notes/iaea/'+fd,
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
                 markersize=3, markeredgewidth=0.75,
                 fillstyle='none', linestyle='none', label='Exp.')
            break
        except:
            pass
                    
    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    title('Target=%s'%a)    
    legend()
    
def gen_plots(zs = dts.ss, mrc=0):
    for z in zs:
        print(z)
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
        prange(z)
        
def pden(z):
    od = 'data/%s'%fac.ATOMICSYMBOL[z]
    d = rden(od)
    clf()
    plot(d[0], d[1])
    if (z == 79):
        r = loadtxt('output/z79s/au_rho.txt', unpack=1)
        plot(r[0], r[1], marker='.')
    xlabel('radius (a.u.)')
    ylabel(r'$4\pi r^2\rho(r)$ (a.u.)')

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
        r = rdedx('d%dB'%i)
        loglog(r[0], r[1], label=r'$\rho$=%d g/cc'%(int(ds[i])))
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')

    legend()
    
    savefig('bd_dedx.png')
    savefig('bd_dedx.pdf')
    
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
        r = rdedx('d%dt%02dB'%(ide,i))
        loglog(r[0], r[1], label=r'T=%g eV'%(int(ts[i])))
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    title(r'$\rho$=%g g/cc'%ds[ide])
    legend()
    
    savefig('bt%d_dedx.png'%(ide))
    savefig('bt%d_dedx.pdf'%(ide))
    
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)
    clf()
    for i in range(nt):
        r = rdedx('d%dt%02dB'%(ide,i))
        loglog(r[0], r[2]/1e3, label=r'T=%g eV'%(int(ts[i])))
    xlabel('Energy (MeV)')         
    ylabel('Proton Range (g/cm$^2$)')          
    title(r'$\rho$=%g g/cc'%ds[ide])
    legend()
    
    savefig('bt%d_range.png'%ide)
    savefig('bt%d_range.pdf'%ide)

def pct(ide=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.5]
    ts = [10, 20, 30]
    nd = len(ds)
    nt = len(ts)
    clf()
    a = 1e-18/(12*1.67e-24*1e6)
    for i in range(nt):
        r = rdedx('d%dt%02dC'%(ide,i))
        plot(r[0], r[1]*a, label=r'T=%g eV'%(int(ts[i])))
    xlim(0, 1.0)
    ylim(0, 1.2)
    d = loadtxt('malko_fig.txt', unpack=1)
    plot(d[0][::2]/1e3, d[1][::2], marker='o', color='k', label='DFT', linestyle='none')
    errorbar(d[0][::2]/1e3, d[1][::2], yerr=(d[1][1::2]-d[1][::2]), capsize=3, marker='o', color='k', linestyle='none')
    xlabel('Energy (MeV)')                   
    ylabel(r'dE/dx (keV/($\mu$g/cm$^2$))')
    title(r'$\rho$=%g g/cc'%ds[0])
    legend()
    
    savefig('ct%d_dedx.png'%(ide))
    savefig('ct%d_dedx.pdf'%(ide))
    
def pnit(ide=0, nde=0):
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
    a = aa.AA()
    for j in range(nd):
        for i in range(nt):
            r = rdedx('d%dt%02dNi'%(j,i))
            fi = interpolate.interp1d(r[0], r[2])
            ry[j,i] = fi(1.6)
            rz[j,i] = a.rden('d%dt%02dNi/Ni'%(j,i), header='zb')
            rf[j,i] = a.rden('d%dt%02dNi/Ni'%(j,i), header='zf')

    for j in range(ide, nde+ide):
        plot(rz[j], ry[j], label=r'$\rho=%g$'%ds[j], marker='o')
        plot(rf[j], ry[j], label=r'$\rho=%g$'%ds[j], marker='o')
    xlabel('Zbar')                   
    ylabel(r'range (mg/cm$^2$)')
    title('Proton in Ni')
    xlim(-1,26)
    ylim(3,20)
    legend()
    
    savefig('nit%d_range.png'%(ide))
    savefig('nit%d_range.pdf'%(ide))
    
def palt(ide=0, nde=0):
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.9,left=0.15,right=0.95)

    ds = [0.027,0.27, 2.7]
    ts = 10**linspace(-1,4.0,35)

    nd = len(ds)
    nt = len(ts)
    if nde == 0:
        nde = nd
    clf()
    ry = zeros((nd,nt))
    rz = zeros((nd,nt))
    rf = zeros((nd,nt))
    a = aa.AA()
    for j in range(nd):
        for i in range(nt):
            r = rdedx('d%dt%02dAl'%(j,i))
            fi = interpolate.interp1d(r[0], r[2])
            ry[j,i] = fi(1.6)
            rz[j,i] = a.rden('d%dt%02dAl/Al'%(j,i), header='zb')
            rf[j,i] = a.rden('d%dt%02dAl/Al'%(j,i), header='zf')

    for j in range(ide,nde+ide):
        plot(rz[j], ry[j], label=r'$\rho=%g$'%ds[j], marker='o')
        plot(rf[j], ry[j], label=r'$\rho=%g$'%ds[j], marker='o')
    xlabel('Zbar')                   
    ylabel(r'range (mg/cm$^2$)')
    title('Proton in Al')
    xlim(-1,12)
    ylim(3,12)
    legend()
    
    savefig('alt%d_range.png'%(ide))
    savefig('alt%d_range.pdf'%(ide))
    
