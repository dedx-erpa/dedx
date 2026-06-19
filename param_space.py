"""
param_space.py -- velocity-ratio / electron-coupling map of reported ion-stopping
experiments, with the 1980s NRL (Young 1982) and Sandia (Olsen 1985) enhanced-
stopping experiments placed on it (the points that were missing from the
compilation).  Editable replacement for the imported map image.

  v_p   = sqrt(2 E / m_proj);  v_th = sqrt(k T_e / m_e)
  Gamma = e^2 / (a_e k T_e),  a_e = (3/4 pi n_e)^(1/3),  n_e = Zbar rho N_A / A

NRL/Sandia plasma conditions (peak-heating case where the enhancement was seen):
  Young 1982: 1-MeV deuterons; Al T_e=13-17 eV, rho~0.02 g/cc; Mylar T_e=9-11 eV.
  Olsen 1985: 1.6-MeV protons; Al T_e~48 eV, Zbar 8-10; Ni T_e~42 eV, Zbar 10-13.
Densities at measurement are expansion-dependent (~1% solid); Gamma ~ rho^(1/3).

The modern literature points are APPROXIMATE, read from the published compilation
(Malko et al. 2022; Frenje et al. 2019) -- replace with exact values / add the
citation before publication.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
# modern literature points (approximate, from the compilation): name, v_p/v_th, Gamma
LIT = [
    ('Frenje 2019', 0.4, 0.011), ('Frenje 2015', 1.4, 0.018), ('Cayzac 2017', 1.6, 0.011),
    ('Frank 2013', 3.2, 0.011), ('Chen 2018', 1.2, 0.006), ('Hicks 2000', 0.8, 0.002),
    ('Couillaud 1994', 2.0, 0.006), ('Sakumi 2001', 3.0, 0.009), ('Jacoby 1995', 2.0, 0.04),
    ('Roth 2000', 4.0, 0.04), ('Flierl 1998', 6.0, 0.03), ('Ren 2020', 8.0, 0.10),
    ('Dietrich 1990', 10.0, 0.10), ('Gardes 1992', 17.0, 0.10), ('Belyaev 1996', 11.0, 0.036),
    ('Hoffmann 1990', 14.0, 0.078), ('Malko 2022', 2.6, 0.70), ('Zylstra 2015', 13.0, 0.50),
    ('Ortner 2015', 18.0, 0.45), ('Lahmann 2023', 33.0, 1.7),
]

fig, ax = plt.subplots(figsize=(9, 6.2))
# ICF alpha-emission band (cold fuel -> hot spot)
from matplotlib.patches import Ellipse
ax.add_patch(Ellipse((0.7, 0.06), width=4.2, height=4.0, angle=33,
                      transform=ax.transData, color='steelblue', alpha=0.12, zorder=0))
ax.scatter([2.0], [0.5], s=170, marker='D', color='blue', zorder=5)
ax.annotate('alpha cold fuel', (2.0, 0.5), (2.4, 0.62), color='blue', fontsize=9, fontweight='bold')
ax.scatter([0.2], [0.01], s=170, marker='D', color='red', zorder=5)
ax.annotate('alpha hot spot', (0.2, 0.01), (0.12, 0.0055), color='red', fontsize=9, fontweight='bold')

# modern literature
for nm, x, y in LIT:
    ax.scatter([x], [y], s=42, marker='o', facecolor='0.75', edgecolor='0.4', zorder=3)
    ax.annotate(nm, (x, y), (x * 1.05, y * 1.12), fontsize=7, color='0.4')

# NRL/Sandia (this paper's addition)
print(' experiment                 | v_p/v_th | Gamma')
for lab, E, m, Te, rho, A, Z in NEW:
    x, g = coords(E, m, Te, rho, A, Z)
    print(' %-26s |  %5.2f   | %.3f' % (lab, x, g))
    ax.scatter([x], [g], s=140, marker='*', color='green', edgecolor='k', zorder=6)
    ax.annotate(lab.split(' (')[0], (x, g), (x * 1.06, g * 1.18), fontsize=8.5,
                color='darkgreen', fontweight='bold')

ax.axvline(1.0, color='0.6', ls=':'); ax.text(1.04, 1.3e-3, 'Bragg peak', fontsize=9)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-1, 1e2); ax.set_ylim(1e-3, 3e0)
ax.set_xlabel(r'velocity ratio  $v_p/v_{th}$', fontsize=12)
ax.set_ylabel(r'electron coupling  $\Gamma$', fontsize=12)
ax.scatter([], [], marker='*', color='green', edgecolor='k', s=140, label='NRL/Sandia (this compilation)')
ax.scatter([], [], marker='o', facecolor='0.75', edgecolor='0.4', s=42, label='reported experiments')
ax.legend(loc='lower right', fontsize=9)
ax.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig('param_space.pdf'); fig.savefig('param_space.png', dpi=130)
print('\nwrote param_space.pdf / .png')
