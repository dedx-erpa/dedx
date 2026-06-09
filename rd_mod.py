#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compound dE/dx comparison plot: overlays PSTAR + IAEA experiment + the
gapless eRPA baseline (data/<a>) + any gap variants (data/<a>_<N>).

Adapted from the Oct 2024 cmp_dedx_mod (dedx_ref) with fixes:
  - qualified rd.e2v / rd.tofloat (were bare names -> NameError)
  - variant blocks collapsed into a loop; labels read each run's mloss
    so gapless (mloss=24, icb=1) and gapped (mloss>=100, icb=0) are
    distinguishable.

@author: mehlhorn
"""
from pylab import *
import matplotlib.pyplot as plt
from pfac import fac
import os
import rd   # need rd functions (rpstar, rdedx, e2v, tofloat)

# suffixes searched for variant runs, in plot order
VARIANTS = ['0', '1', '2', '3', '4', '23', '101']


def _label(od, fallback):
    """Label a run by its loss mode; flag gap (icb=0) runs."""
    try:
        h = rd.rdedx(od, header='')
        ml = int(h['mloss'])
        if ml >= 100:
            return 'gap (mloss=%d)' % ml
        return 'eRPA (mloss=%d)' % ml
    except Exception:
        return fallback


def cmp_dedx_mod(zt, xsc='e', pstar=1, exp=1, orpa=1):
    clf()
    plt.rcParams.update({'font.size': 15})
    plt.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.95)
    if type(zt) == type(0):
        a = fac.ATOMICSYMBOL[zt]
    else:
        a = zt
        try:
            zt = fac.ATOMICSYMBOL.index(a)
        except Exception:
            zt = 0
    print('a = ', a)
    d = 'data/%s' % a

    if xsc == 'e':
        xlabel('Energy (MeV)')
        xlim(1e-3, 1e2)
    else:
        xlabel('Velocity (a.u.)')
        xlim(0, 20.0)

    def _plot(x, y, **kw):
        if xsc == 'v':
            plot(rd.e2v(x), y, **kw)
        else:
            semilogx(x, y, **kw)

    # --- PSTAR ---
    if pstar and os.path.exists(d + '/pstar.txt'):
        r = rd.rpstar(d)
        _plot(r[0], r[1], marker='o', markerfacecolor='none',
              color='k', label='PSTAR')

    # --- IAEA experiment ---
    if exp:
        fds = [a + '.txt', a.lower() + '.txt', a.upper() + '.txt',
               'h' + a.lower() + '.txt', 'H' + a + '.txt', 'H' + a.upper() + '.txt']
        if a == 'H4C5O2':
            fds.append('MYLAR.txt')
        for fd in fds:
            try:
                r = loadtxt('data/iaea/' + fd, unpack=1,
                            usecols=(0, 1), dtype=str)
                x = array([rd.tofloat(r[0][i]) for i in range(len(r[0]))])
                y = array([rd.tofloat(r[1][i]) for i in range(len(r[1]))])
                if ((zt != 14 and x[0] == 0 and y[0] == 0 and x[1] == 0) or
                        max(y) < 0.1):
                    y = y * fac.ATOMICMASS[zt] * 1.67
                else:
                    x = x / 1e3
                plot(x, y, marker='o', color='r', markersize=3,
                     markeredgewidth=0.5, fillstyle='none',
                     linestyle='none', label='Exp.')
                break
            except Exception:
                pass

    # --- gapless baseline: data/<a> ---
    if os.path.exists(d + '/dedx.dat'):
        y = rd.rdedx(d)
        _plot(y[0], y[1], linewidth=2, label=_label(d, 'eRPA'))

    # --- gap variants: data/<a>_<N> ---
    for n in VARIANTS:
        od = '%s_%s' % (d, n)
        if os.path.exists(od + '/dedx.dat'):
            y = rd.rdedx(od)
            _plot(y[0], y[1], linestyle='--', linewidth=3,
                  label=_label(od, 'variant_%s' % n))

    # --- Wang et al overlay (optional) ---
    if orpa and os.path.exists(d + '/odedx.dat'):
        y = loadtxt(d + '/odedx.dat', unpack=1)
        x = y[0]
        if max(x) > 1000:
            x /= 1e3
        plot(x, y[1], linestyle='dashdot', linewidth=3, label='Wang et al')

    ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    title('Target=%s' % a)
    legend()
