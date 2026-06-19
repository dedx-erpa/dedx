"""
param_space.py -- velocity-ratio / electron-coupling map of reported ion-stopping
experiments, recreating Fig. 1 of the GSI flagship proposal "Experiments on energy
loss of ions in fusion relevant extreme states of matter" (Volpe et al., 2024) and
OVERLAYING the 1980s NRL (Young 1982) and Sandia (Olsen 1985) enhanced-stopping
experiments -- the first measurements of enhanced stopping in partially ionized
plasmas, which the modern compilation (and the proposal's Fig. 1) omit.

  v_p   = sqrt(2 E / m_proj);  v_th = sqrt(k T_e / m_e)
  Gamma = e^2 / (a_e k T_e),  a_e = (3/4 pi n_e)^(1/3),  n_e = Zbar rho N_A / A

This uses the SAME definitions as the proposal (its Fig. 3 box):
  Gamma = q^2/(r_s k_B T),  r_s = (3/4 pi n)^(1/3).

NRL/Sandia plasma conditions (peak-heating case where the enhancement was seen):
  Young 1982: 1-MeV deuterons; Al T_e=13-17 eV, rho~0.02 g/cc; Mylar T_e=9-11 eV.
  Olsen 1985: 1.6-MeV protons; Al T_e~48 eV, Zbar 8-10; Ni T_e~42 eV, Zbar 10-13.
Densities at measurement are expansion-dependent (~1% solid); Gamma ~ rho^(1/3).
Computed coords: Young Al (6.0, 0.17), Young Mylar (7.4, 0.30),
                 Olsen Al (6.0, 0.085), Olsen Ni (6.4, 0.12)
-> all four land INSIDE the proposal's "cold fuel" / unexplored target box.

Literature points and their plasma-generation method are read from the proposal's
Fig. 1 (Volpe et al. 2024), which is itself an update of the compilation in
Malko et al., Nat. Commun. 13, 2893 (2022), Fig. 1 (+ Frenje et al. 2019).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

c, mec2, e2, NA = 2.998e10, 511000.0, 1.44e-7, 6.022e23

def coords(E_MeV, m_amu, Te, rho, A, Z):
    vp = c * np.sqrt(2 * E_MeV / (m_amu * 931.494))
    vth = c * np.sqrt(Te / mec2)
    ne = Z * rho * NA / A
    ae = (3 / (4 * np.pi * ne)) ** (1 / 3.)
    return vp / vth, e2 / (ae * Te)

# this work: NRL/Sandia  (label, E, m, Te, rho, A, Z)
NEW = [
    ('Young 1982 Al (NRL)',    1.0, 2.014, 15, 0.020, 27.0, 3.0),
    ('Young 1982 Mylar (NRL)', 1.0, 2.014, 10, 0.010, 8.73, 3.0),
    ('Olsen 1985 Al (Sandia)', 1.6, 1.007, 48, 0.027, 27.0, 9.0),
    ('Olsen 1985 Ni (Sandia)', 1.6, 1.007, 42, 0.089, 58.7, 11.5),
]

# plasma-generation method -> color/label (matches the proposal legend)
METHODS = {
    'laser':    ('#1f6fb4', 'laser plasma'),
    'xray':     ('#8aa829', 'X-ray driven'),
    'zpinch':   ('#e8a200', 'Z-pinch'),
    'gas':      ('#b81d2a', 'gas discharge'),
    'exploding':('#6a2c8f', 'exploding pusher'),
}

# modern literature points: name, v_p/v_th, Gamma, method, (label dx, dy mult)
# positions digitized from the proposal's Fig. 1 (Volpe et al. 2024)
LIT = [
    ('Lahmann 2023',  33.0, 2.0,  'xray',      (0.55, 1.0)),
    ('Zylstra 2015',  20.0, 1.0,  'xray',      (1.05, 1.0)),
    ('Ortner 2015',   15.0, 0.5,  'xray',      (1.05, 1.0)),
    ('Malko 2022',     4.5, 0.7,  'laser',     (0.30, 1.1)),
    ('Ren 2020',       8.0, 0.11, 'xray',      (0.45, 1.25)),
    ('Dietrich 1990', 10.0, 0.11, 'zpinch',    (1.07, 1.2)),
    ('Hoffmann 1990', 11.0, 0.07, 'zpinch',    (0.40, 1.25)),
    ('Gardes 1992',   18.0, 0.085,'gas',       (1.07, 1.0)),
    ('Belyaev 1996',  16.0, 0.06, 'gas',       (1.07, 0.85)),
    ('Flierl 1998',   11.0, 0.04, 'zpinch',    (1.07, 1.0)),
    ('Roth 2000',      5.0, 0.04, 'laser',     (1.07, 1.0)),
    ('Jacoby 1995',    3.0, 0.05, 'gas',       (0.35, 1.1)),
    ('Frenje 2015',    1.6, 0.02, 'exploding', (0.40, 1.0)),
    ('Cayzac 2017',    1.3, 0.011,'laser',     (0.42, 1.15)),
    ('Sakumi 2001',    2.5, 0.011,'laser',     (1.07, 1.0)),
    ('Frank 2013',     3.0, 0.007,'laser',     (1.07, 1.0)),
    ('Couillaud 1994', 2.0, 0.005,'laser',     (1.07, 0.85)),
    ('Frenje 2019',    1.0, 0.01, 'exploding', (0.30, 1.0)),
    ('Chen 2018',      1.2, 0.0025,'laser',    (1.07, 1.0)),
    ('Hicks 2000',     0.8, 0.002,'exploding', (0.30, 1.0)),
]

fig, ax = plt.subplots(figsize=(9.2, 6.4))

def log_ellipse(ax, xc, yc, wdex, hdex, angle_deg, **kw):
    """Ellipse drawn in log10 space so it stays elliptical on log-log axes.
    wdex/hdex = full axis lengths in decades; angle in log-space degrees."""
    t = np.linspace(0, 2 * np.pi, 120)
    a, b = wdex / 2., hdex / 2.
    th = np.radians(angle_deg)
    lx = np.log10(xc) + a * np.cos(t) * np.cos(th) - b * np.sin(t) * np.sin(th)
    ly = np.log10(yc) + a * np.cos(t) * np.sin(th) + b * np.sin(t) * np.cos(th)
    ax.fill(10 ** lx, 10 ** ly, **kw)

# --- proposal regions ----------------------------------------------------
# alpha hot spot: orange ellipse, lower-left (drawn in log space)
log_ellipse(ax, 0.30, 0.045, 0.95, 0.75, 33, color='#f0a030', alpha=0.22, zorder=0)
ax.text(0.20, 0.05, r'$\alpha$ hot spot', color='#c06000', fontsize=11,
        fontstyle='italic', ha='center', va='center', zorder=1)
# cold fuel / proposal's unexplored WDM target box: light-blue rectangle
ax.add_patch(Rectangle((1.0, 0.1), 9.0, 1.9, transform=ax.transData,
                       color='#9cc3e0', alpha=0.28, zorder=0))
ax.text(2.3, 1.4, r'$\alpha$ cold fuel', color='#2a6699', fontsize=11,
        fontstyle='italic', ha='center', va='center', zorder=1)
ax.text(3.0, 0.22, 'unexplored WDM\n(this proposal)', color='#2a6699', fontsize=8.5,
        ha='center', va='center', alpha=0.9, zorder=1)

# --- literature, colored by plasma-generation method ---------------------
for nm, x, y, meth, (dx, dy) in LIT:
    col = METHODS[meth][0]
    ax.scatter([x], [y], s=70, marker='o', color=col, edgecolor='0.25',
               linewidth=0.4, zorder=4)
    ax.annotate(nm, (x, y), (x * dx, y * dy), fontsize=7.5, color='0.15')

# --- NRL/Sandia (this work's addition): green stars ----------------------
# hand-placed label offsets (xytext in data coords) so the cluster reads cleanly
NEW_LBL = {
    'Young 1982 Al':    (3.7, 0.24),
    'Young 1982 Mylar': (8.0, 0.40),
    'Olsen 1985 Al':    (3.0, 0.072),
    'Olsen 1985 Ni':    (7.2, 0.135),
}
print(' experiment                 | v_p/v_th | Gamma')
for lab, E, m, Te, rho, A, Z in NEW:
    x, g = coords(E, m, Te, rho, A, Z)
    print(' %-26s |  %5.2f   | %.3f' % (lab, x, g))
    ax.scatter([x], [g], s=210, marker='*', color='#2ca02c', edgecolor='k',
               linewidth=0.8, zorder=7)
    short = lab.split(' (')[0]
    tx, ty = NEW_LBL[short]
    ax.annotate(short, (x, g), (tx, ty), fontsize=8.5, color='#1a661a',
                fontweight='bold', zorder=7,
                arrowprops=dict(arrowstyle='-', color='#1a661a', lw=0.5))

# --- axes / Bragg peak ---------------------------------------------------
ax.axvline(1.0, color='0.6', ls=':'); ax.text(1.04, 1.4e-3, 'Bragg peak', fontsize=9)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-1, 1e2); ax.set_ylim(1e-3, 3e0)
ax.set_xlabel(r'velocity ratio  $v_p/v_{th}$', fontsize=12)
ax.set_ylabel(r'electron coupling  $\Gamma$', fontsize=12)

# --- legends -------------------------------------------------------------
meth_handles = [plt.Line2D([], [], marker='o', ls='', color=col, mec='0.25',
                           ms=8, label=lab) for col, lab in METHODS.values()]
nrl_handle = plt.Line2D([], [], marker='*', ls='', color='#2ca02c', mec='k',
                        ms=14, label='NRL / Sandia (this work)')
leg1 = ax.legend(handles=meth_handles, loc='upper left', fontsize=8.5,
                 framealpha=0.9, title='plasma generation')
ax.add_artist(leg1)
ax.legend(handles=[nrl_handle], loc='lower right', fontsize=9, framealpha=0.9)

ax.grid(alpha=.3, which='both')
fig.tight_layout()
fig.savefig('param_space.pdf'); fig.savefig('param_space.png', dpi=140)
print('\nwrote param_space.pdf / .png')
