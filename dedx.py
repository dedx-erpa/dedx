from pfac import aa
from pfac import fac
import rd
from optparse import OptionParser
import os, time
from dts import *

def prep_inp(opts, od, zt):
    with open('%s/odir.inp'%od, 'w') as f:
        f.write('"%s"\n'%od)
        if opts.floss == '':
            f.write('NONE\n')
        else:
            f.write('"%s"\n'%opts.floss)
        if opts.floss0 == '':
            f.write('NONE\n')
        else:
            f.write('"%s"\n'%opts.floss0)

    with open('%s/dedx.inp'%od, 'w') as f:
        f.write('&dedxinp\n\n')
        f.write('  zzp = %g\n'%opts.zp)
        f.write('qmass = %9.3e\n'%fac.ATOMICMASS[int(opts.zp)])
        f.write('  ztg = %d\n'%zt)
        f.write('amass = %9.3e\n'%fac.ATOMICMASS[zt])
        f.write('  mep = %d\n'%opts.mep)
        f.write('  emin = %9.3e\n'%opts.emin)
        f.write('  emax = %9.3e\n'%opts.emax)
        f.write('  mloss = %d\n'%opts.mloss)
        f.write('  mout = %d\n'%opts.mout)
        f.write('  dinp = %9.3e\n'%opts.d)
        f.write('  tinp = %9.3e\n'%opts.t)
        f.write('  epa = %9.3e\n'%opts.epa)
        f.write('  epb = %9.3e\n'%opts.epb)
        f.write('  epc = %9.3e\n'%opts.epc)
        f.write('  epd = %9.3e\n'%opts.epd)
        f.write('  epe = %9.3e\n'%opts.epe)
        f.write('&end\n')

    if (opts.frho != '' and
        opts.frho != 'NONE' and
        os.path.exists(opts.frho)):
        os.system('cp %s %s/rho.functions'%(opts.frho, od))
        opts.aa = 0
        
p = OptionParser()
p.add_option('--zp', dest='zp', type='float',
             default=1, help='projectile z')
p.add_option('--zt', dest='zt', type='int',
             default=13, help='target z')
p.add_option('--zc', dest='zc', type='string',
             default='', help='compound z array')
p.add_option('--wc', dest='wc', type='string',
             default='', help='compound w array')
p.add_option('--fc', dest='fc', type='string',
             default='', help='compound formula')
p.add_option('--d', dest='d', type='float',
             default=0.0, help='density in g/cc')
p.add_option('--t', dest='t', type='float',
             default=1.0, help='temperature in eV')
p.add_option('--taa', dest='taa', type='float',
             default=0.5, help='min temp in eV to run avg atom model')
p.add_option('--frho', dest='frho', type='string',
             default='', help='rho file to copy')
p.add_option('--od', dest='od', type='string',
             default='data', help='output directory')
p.add_option('--emin', dest='emin', type='float',
             default=1e-3, help='min energy in MeV')
p.add_option('--emax', dest='emax', type='float',
             default=100.0, help='max energy in MeV')
p.add_option('--mep', dest='mep', type='int',
             default=100, help='number of energy points')
p.add_option('--mloss', dest='mloss', type='int',
             default=24, help='energy loss func mode')
p.add_option('--floss', dest='floss', type='string',
             default='NONE', help='energy loss func data file')
p.add_option('--floss0', dest='floss0', type='string',
             default='NONE', help='2nd energy loss func data file')
p.add_option('--aa', dest='aa', type='int',
             default=2, help='generate aa density distribution')
p.add_option('--bqp', dest='bqp', type='float',
             default=-1e12, help='aa boundary condition')
p.add_option('--nr', dest='nr', type='int',
             default=600, help='number of radial grid')
p.add_option('--rmin', dest='rmin', type='float',
             default=7.5e-5, help='minimum radial grid')
p.add_option('--dedx', dest='dedx', type='string',
             default='./dedx', help='dedx executable')
p.add_option('--mout', dest='mout', type='int',
             default=0, help='output radial dist. of dedx')
p.add_option('--sc', dest='sc', type='int',
             default=0, help='sc param of aa.AA')
p.add_option('--pmi', dest='pmi', type='int',
             default=0, help='pmi param of aa.AA')
p.add_option('--v', dest='v', type='int',
             default=0, help='verbose level')
p.add_option('--maa', dest='maa', type='int',
             default=1, help='mode of compound aa model')
p.add_option('--npaa', dest='npaa', type='int',
             default=0, help='num. proc for aa model')
p.add_option('--epa', dest='epa', type='float',
             default=-1e11, help='epa param')
p.add_option('--epb', dest='epb', type='float',
             default=-1e11, help='epb param')
p.add_option('--epc', dest='epc', type='float',
             default=-1e11, help='epc param')
p.add_option('--epd', dest='epd', type='float',
             default=-1e11, help='epd param')
p.add_option('--epe', dest='epe', type='float',
             default=-1e11, help='epe param')

opts,args = p.parse_args()

if not os.path.exists(opts.od):
    os.system('mkdir -p %s'%opts.od)
    
fn = '%s/dedx.par'%opts.od
with open(fn, 'w') as f:
    for kw in opts.__dict__.keys():
        f.write('%s = %s\n'%(kw, str(opts.__dict__[kw])))
        
ts = 10**numpy.linspace(-1.,4.,21)
if (opts.floss == 'NONE' and
    opts.floss0 == 'NONE' and
    abs(opts.mloss)%10 != 0):
    w = numpy.where(ts > opts.t)[0]
    if len(w) == 0:
        opts.floss = 't20.dat'
    else:
        w = w[0]
        opts.floss = 't%02d.dat'%w
        if w > 0:
            opts.floss0 = 't%02d.dat'%(w-1)
        else:
            opts.floss0 = 'tct.dat'

if opts.zc != '':
    zc = [int(x) for x in opts.zc.split(',')]
    wc = [float(x) for x in opts.wc.split(',')]
elif opts.fc != '':
    zc,wc = aa.zw4c(opts.fc)
    if (len(zc) == 1):
        zt = zc[0]
        zc = []
        wc = []
else:
    zc = []
    wc = []    

if (opts.d <= 0):
    if opts.fc != '':
        opts.d = getden(opts.fc)
    elif opts.zt > 0:
        opts.d = getden(opts.zt)
if opts.aa > 1:
    t0 = time.time()
    taa = max(opts.t, opts.taa)
    if len(zc) > 0:
        if opts.v > 0:
            print('running avgatom zc=%s wc=%s fc=%s d=%g t=%g'%(opts.zc, opts.wc, opts.fc, opts.d, taa))
        a = aa.AA(z=zc, d=opts.d, t=taa, wm=wc, dd=opts.od, bqp=opts.bqp, sc=opts.sc, pmi=opts.pmi, nc=opts.npaa)
        a.run(opts.maa)
        t1 = time.time()
        if opts.v > 0:
            print('done avgatom zc=%s wc=%s fc=%s d=%g t=%g in %10.3E'%(opts.zc, opts.wc, opts.fc, opts.d, taa,(t1-t0)))
    else:
        if opts.v > 0:
            print('running avgatom z=%d d=%g t=%g ...'%(opts.zt,opts.d,taa))
        a = aa.AA(z=opts.zt, d=opts.d, t=taa, dd=opts.od, bqp=opts.bqp, sc=opts.sc, pmi=opts.pmi)
        a.run()
        t1 = time.time()
        if opts.v > 0:
            print('done avgatom z=%d d=%g t=%g in %10.3E'%(opts.zt, opts.d, taa,(t1-t0)))
    
if opts.aa > 0:
    if opts.aa == 1:
        a = aa.AA()
    if len(zc) > 0:
        for z in zc:
            az = fac.ATOMICSYMBOL[z]
            od = '%s/%s'%(opts.od,az)
            if not os.path.exists(od):
                os.system('mkdir -p %s'%od)
            a.wden(od, opts.nr, '%s/rho.functions'%od,
                   rmin=opts.rmin/opts.zt**2)
    else:
        a.wden('%s/%s'%(opts.od,fac.ATOMICSYMBOL[opts.zt]),
               opts.nr,
               '%s/rho.functions'%opts.od,
               rmin=opts.rmin/opts.zt**2)
        if opts.v > 0:
            print('rho.functions written')

if len(zc) > 0:
    for z in zc:
        az = fac.ATOMICSYMBOL[z]
        od = '%s/%s'%(opts.od,az)
        prep_inp(opts, od, z)
else:
    prep_inp(opts, opts.od, opts.zt)

if opts.dedx != '' and opts.dedx != 'None':
    t0 = time.time()
    if opts.v > 0:
        print('running dedx ...')
    if len(zc) > 0:
        rx = None
        ry = None
        wt = 0.0
        am = 0.0
        azt = 0.0
        azb = 0.0
        ars = 0.0
        for iz in range(len(zc)):
            z = zc[iz]
            w = wc[iz]
            az = fac.ATOMICSYMBOL[z]
            am = am + fac.ATOMICMASS[z]*w
            azt = azt + z*w
            od = '%s/%s'%(opts.od,az)
            os.system('%s < %s/odir.inp'%(opts.dedx, od))
            h = rd.rdedx(od, header='')
            r = rd.rdedx(od)
            zb = h['zbar']
            azb = azb + zb*w
            ars = ars + h['rs']*w
            if rx is None:
                rx = r[0]
                ry = r[1]*w
                wt = w
            else:
                ry += r[1]*w
                wt += w
        ry /= wt
        am /= wt
        azt /= wt
        azb /= wt
        ars /= wt
        rz = rd.int_range(numpy.array((rx,ry)),m=am)
        with open('%s/dedx.dat'%(opts.od), 'w') as f:
            f.write('#   nzt = %d\n'%len(zc))
            sz = '#    zt ='
            sw = '#    wt ='            
            for iz in range(len(zc)):
                sz += ' %10.4E'%zc[iz]
                sw += ' %10.4E'%wc[iz]
            f.write(sz+'\n')
            f.write(sw+'\n')
            f.write('#    zp = %10.4E\n'%opts.zp)
            f.write('#    rs = %15.8E\n'%ars)
            f.write('#    Te = %15.8E\n'%opts.t)
            f.write('#   rho = %15.8E\n'%opts.d)
            f.write('#  zbar = %15.8E\n'%azb)
            f.write('#   mep = %3d\n'%len(rx))
            f.write('# mloss = %3d\n'%opts.mloss)
            f.write('#   E/AMU (MeV)            dEdX           range\n')
            for i in range(len(rx)):
                f.write('%15.8E %15.8E %15.8E\n'%(rx[i],ry[i],rz[i]))
    else:
        os.system('%s < %s/odir.inp'%(opts.dedx,opts.od))
    t1 = time.time()
    if opts.v > 0:
        print('done in %10.3E s'%(t1-t0))
