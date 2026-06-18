"""
param_space.py -- place the 1980s NRL (Young 1982) and Sandia (Olsen 1985)
enhanced-stopping experiments on the velocity-ratio / electron-coupling map
(v_p/v_th vs Gamma), the parameter space used by Frenje 2019 / Malko 2022.

Coordinates:
  v_p   = sqrt(2 E / m_proj)                 (projectile velocity)
  v_th  = sqrt(k T_e / m_e)                  (electron thermal velocity)
  Gamma = e^2 / (a_e k T_e),  a_e = (3/4 pi n_e)^(1/3),  n_e = Zbar rho N_A / A

Plasma conditions (from the papers, peak heating / high-intensity case where the
enhancement was measured):
  Young 1982 (NRL, 1-MeV deuterons, 250 kA/cm^2): Al  T_e=13-17 eV, rho~0.02 g/cc;
                                                  Mylar T_e=9-11 eV, rho~0.01 g/cc.
  Olsen 1985 (Sandia, 1.6-MeV protons, 1.4 TW/cm^2): Al T_e~48 eV, Zbar 8-10;
                                                     Ni T_e~42 eV, Zbar 10-13.
Densities at the measurement are expansion-dependent; we take ~1% solid (consistent
with Young's stated 0.02 g/cc Al).  Gamma ~ rho^(1/3), so a 10x density change moves
Gamma by ~2.2x; v_p/v_th is density-independent.
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

# (label, E, m_proj, Te, rho, A_target, Zbar, marker)
EXPTS = [
    ('Young 1982 Al (NRL)',    1.0, 2.014, 15, 0.020, 27.0, 3.0, '*'),
    ('Young 1982 Mylar (NRL)', 1.0, 2.014, 10, 0.010, 8.73, 3.0, '*'),
    ('Olsen 1985 Al (Sandia)', 1.6, 1.007, 48, 0.027, 27.0, 9.0, 'D'),
    ('Olsen 1985 Ni (Sandia)', 1.6, 1.007, 42, 0.089, 58.7, 11.5, 'D'),
]

fig, ax = plt.subplots(figsize=(8.5, 6))
# ICF alpha reference zones (read from the Frenje/Malko map)
ax.scatter([2.0], [0.5], s=160, marker='D', color='blue', zorder=5)
ax.annotate('alpha cold fuel', (2.0, 0.5), (2.3, 0.6), color='blue', fontsize=9)
ax.scatter([0.2], [0.01], s=160, marker='D', color='red', zorder=5)
ax.annotate('alpha hot spot', (0.2, 0.01), (0.12, 0.006), color='red', fontsize=9)

print(' experiment                 | v_p/v_th | Gamma')
for lab, E, m, Te, rho, A, Z, mk in EXPTS:
    x, g = coords(E, m, Te, rho, A, Z)
    print(' %-26s |  %5.2f   | %.3f' % (lab, x, g))
    ax.scatter([x], [g], s=110, marker=mk, color='green', edgecolor='k', zorder=6)
    ax.annotate(lab.split(' (')[0], (x, g), (x * 1.05, g * 1.12), fontsize=8.5, color='darkgreen')

ax.axvline(1.0, color='0.7', ls=':'); ax.text(1.02, 1.2e-3, 'Bragg peak', fontsize=9)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-1, 1e2); ax.set_ylim(1e-3, 1e0)
ax.set_xlabel(r'velocity ratio  $v_p/v_{th}$'); ax.set_ylabel(r'electron coupling  $\Gamma$')
ax.set_title('1980s NRL/Sandia enhanced-stopping experiments in the\n'
             'velocity-ratio / electron-coupling parameter space')
ax.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig('param_space.pdf'); fig.savefig('param_space.png', dpi=130)
print('\nwrote param_space.pdf / .png')
