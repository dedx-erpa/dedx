"""
exp_zylstra.py -- compare the total eRPA model to Zylstra 2015 (PRL 114, 215002):
energy downshift (-Delta E) of 14.7 MeV D3He protons through a Be foil, cold vs
isochorically-heated warm-dense (Te~32 eV).  Zylstra Fig 4 data digitized in
data/refs/Zylstra_fig4_{cold,warm}.csv.

Be foil: rho = 1.77 g/cc, thickness 532 um -> areal density 94 mg/cm^2
(6.29e21 atoms/cm^2) -- a flat foil, so this is a PARAMETER-FREE absolute test.

Model dirs:
  python dedx.py --zt=4 --zp=1 --d=1.85 --t=32    --aa=2 --nuc=1 --npot=gk --od=/tmp/zylstra_wdm  --emin=0.3 --emax=15 --mep=45
  python dedx.py --zt=4 --zp=1 --d=1.85 --t=0.025 --aa=2 --nuc=1 --npot=gk --od=/tmp/zylstra_cold --emin=0.3 --emax=15 --mep=45
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

E0 = 14.7                      # D3He proton birth energy (MeV)
NA = 6.022e23
A_Be = 9.012
N_areal = (1.77 / A_Be) * NA * 532e-4    # atoms/cm^2  (= 6.29e21)
print('Be areal density = %.3e atoms/cm^2 (%.0f mg/cm^2)'
      % (N_areal, 1.77 * 532e-4 * 1e3))


def dEloss(od, E0=E0, nstep=400):
    """integrate proton slowing-down through the foil -> -Delta E [MeV]."""
    d = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
    E, dt = d[:, 0], d[:, 3]                # E/AMU [MeV], total stopping [1e-15 eV cm2/atom]
    loge = np.log(E)
    eps = lambda EA: np.exp(np.interp(np.log(max(EA, E[0])), loge, np.log(dt)))
    Emev = E0
    dN = N_areal / nstep
    for _ in range(nstep):
        Emev = max(Emev - eps(Emev) * 1e-15 * dN * 1e-6, 1e-3)
    return E0 - Emev


# model predictions (absolute)
dE_cold = dEloss('/tmp/zylstra_cold')
dE_warm = dEloss('/tmp/zylstra_wdm')
print('model -dE: cold = %.3f MeV, warm = %.3f MeV (warm/cold = %.3f)'
      % (dE_cold, dE_warm, dE_warm / dE_cold))

# digitized Zylstra Fig 4 (x category index, dE MeV; rows in groups of 1 or 3)
cold = np.loadtxt('data/refs/Zylstra_fig4_cold.csv', delimiter=',')
warm = np.loadtxt('data/refs/Zylstra_fig4_warm.csv', delimiter=',')


def group(arr):
    out = []
    for xc in sorted(set(np.round(arr[:, 0]))):
        m = np.abs(arr[:, 0] - xc) < 0.5
        ys = arr[m, 1]
        out.append((np.median(arr[m, 0]), np.median(ys), ys.min(), ys.max()))
    return np.array(out)


cg, wg = group(cold), group(warm)
# category labels (from the figure)
cold_lab = ['72025*', '72026*', 'SRIM', 'ICRU', 'AA-LDA']
warm_lab = ['72018*', '72024*', 'B+F', 'AA-LDA', 'CIP']

fig, (axc, axw) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, g, labs, dEm, title in [(axc, cg, cold_lab, dE_cold, 'Cold Be'),
                                (axw, wg, warm_lab, dE_warm, 'Warm Be (Te≈32 eV)')]:
    x = np.arange(len(g))
    isdata = np.array(['*' in l for l in labs])
    err = [g[:, 1] - g[:, 2], g[:, 3] - g[:, 1]]
    ax.errorbar(x[isdata], g[isdata, 1], yerr=[err[0][isdata], err[1][isdata]],
                fmt='ks', ms=8, capsize=4, label='data (shots)')
    ax.errorbar(x[~isdata], g[~isdata, 1], yerr=[err[0][~isdata], err[1][~isdata]],
                fmt='o', mfc='none', mec='0.4', ms=8, capsize=4, label='other theories')
    ax.axhline(dEm, color='r', lw=2, label='this work (eRPA-AA+nuc) = %.2f' % dEm)
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=35, fontsize=8)
    ax.set_title(title); ax.grid(alpha=.3, axis='y')
axc.set_ylabel(r'$-\Delta E$ (MeV)')
axc.legend(fontsize=8, loc='lower right')
fig.suptitle('Zylstra 2015: 14.7 MeV proton downshift through Be (94 mg/cm², parameter-free)')
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('exp_zylstra.pdf'); fig.savefig('exp_zylstra.png', dpi=130)
print('wrote exp_zylstra.pdf / .png')
print('\n cold: data~%.2f, AA-LDA(Zylstra)~%.2f, this work=%.2f'
      % (cg[:2, 1].mean(), cg[-1, 1], dE_cold))
print(' warm: data~%.2f, AA-LDA(Zylstra)~%.2f, this work=%.2f'
      % (wg[:2, 1].mean(), wg[3, 1], dE_warm))
