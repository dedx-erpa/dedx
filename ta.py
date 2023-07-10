import dief
from pylab import *
import matplotlib.pyplot as plt

def getdata():
    d0 = 1e25
    t0 = 1e3
    es = linspace(0.1, 20.0, 20)
    ep = es/4.0*1e6
    vp = sqrt(ep*2/27.2/1823.)
    r0 = dief.dedxmp(ep, d0, t0, mgk=1, zp=2.0, np = 20)
    a = 4*4*pi/vp**2*27.2*1e-6*0.53e-8**2*d0*1e-4
    y0 = (r0[1]+r0[3])*a
    y1 = y0 + r0[2]*a

    e0 = 3.5
    ep0 = e0/4.0*1e6
    vp0 = sqrt(ep0*2/27.2/1823.)
    a0 = 4*4*pi/vp0**2*27.2*1e-6*0.53e-8**2*d0*1e-4
    ts = 10**(linspace(0.0, 3.0, 20))
    ds = 10**(linspace(24.0, 28.0, 20))
    rt = dief.dedxmp(ep0, d0, ts, mgk=1, zp=2.0, np=20)
    yt0 = (rt[1]+rt[3])*a0
    yt1 = yt0 + rt[2]*a0

    rd = dief.dedxmp(ep0, ds, t0, mgk=1, zp=2.0, np=20)
    yd0 = (rd[1]+rd[3])*a0
    yd1 = yd0 + rd[2]*a0

    zs = linspace(2., 40, 20)
    rz = dief.dedxmp(ep0, d0, t0, mgk=1, zp=zs, np=20)
    yz0 = (rz[1]+rz[3])*a0/4
    yz1 = yz0 + rz[2]*a0/4

    return (es,y0,y1),(ts,yt0,yt1),(ds,yd0,yd1),(zs,yz0,yz1)

def pdata(r, m):
    clf()
    plt.rcParams.update({'font.size':15})
    plt.subplots_adjust(bottom=0.15,top=0.95,left=0.15,right=0.95)
    d = loadtxt('wf%d.txt'%m, unpack=1)
    plot(r[m-1][0], r[m-1][1], label='eRPA-LDA')
    plot(d[0], d[1], marker='o', linestyle=' ', label='qGD')
    legend()
    if m == 2 or m == 3:
        yscale('log')
        xscale('log')
    if (m == 1):
        xlim(-0.5,20.5)
        xlabel(r'$E_p$ (MeV)')
        ylabel(r'dE/dx (Mev/$\mu$m)')
        savefig('dedx_ep.png')
    elif (m == 2):
        xlim(5,1e3)
        ylim(0.5,10)
        xlabel(r'$T_e$ (eV)')
        ylabel(r'dE/dx (Mev/$\mu$m)')
        savefig('dedx_te.png')
    elif (m == 3):
        xlim(1e25,1e28)
        xlabel(r'$n_e$ (cm$^{-3}$)')
        ylabel(r'(n0/n)dE/dx (Mev/$\mu$m)')
        savefig('dedx_ne.png')

if __name__ == '__main__':
    r = getdata()
    for m in [1,2,3]:
        pdata(r, m)
        
