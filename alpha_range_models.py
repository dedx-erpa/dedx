# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""alpha_range_models.py -- DT-alpha range vs temperature: eRPA (this work) vs the
published MD/BPS/Atzeni parameterizations, showing the low-velocity divergence.
eRPA points computed with dedx.py at ne-matched H density; fits from Zylstra &
Hurricane PoP 26,062701(2019) Eq.12/Tab.I (MD,BPS), Atzeni Eq.8. BPS-theory point
= digamma-corrected Coulomb log (Born regime). 3.5 MeV alpha in DT, 10 g/cc."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
c=2.998e10; mec2=511000.0
def MD(T,r): rr=r/100.; return (0.01705*T**1.0706/(1+0.01133*T**1.0706))*(1+1.5793*rr**0.1987)
def BPSf(T,r): rr=r/100.; return (0.01556*T**1.0804/(1+0.007144*T**1.0804))*(1+1.7823*rr**0.2003)
def Atz(T): return 0.025*T**1.25/(1+0.0082*T**1.25)
rho=10.0
Tg=np.linspace(1,10,60)
# eRPA computed points (this run)
Te=np.array([1,2,5,10]); eR=np.array([0.035,0.078,0.244,0.599])
vp=c*np.sqrt(2*0.875/931.494); vr=vp/(c*np.sqrt(Te*1e3/mec2))
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12,5))
ax.plot(Tg,[MD(t,rho) for t in Tg],'b-',label='Zylstra MD fit')
ax.plot(Tg,[BPSf(t,rho) for t in Tg],'b--',label='Zylstra BPS fit')
ax.plot(Tg,[Atz(t) for t in Tg],'g-.',label='Atzeni')
ax.plot(Te,eR,'ro-',ms=8,label='eRPA (this work)')
ax.plot([10],[0.69],'k*',ms=16,label='BPS theory (digamma)')
ax.set_xlabel('T (keV)'); ax.set_ylabel(r'alpha range $\rho\lambda$ (g/cm$^2$)')
ax.set_title('DT 10 g/cc: alpha range vs T'); ax.legend(fontsize=9); ax.grid(alpha=.3)
# ratio vs velocity ratio
ax2.plot(vr,eR/np.array([MD(t,rho) for t in Te]),'ro-',ms=8)
ax2.axhline(1,color='0.5',ls=':'); ax2.axvline(1,color='0.7',ls=':')
for x,t in zip(vr,Te): ax2.annotate('%d keV'%t,(x,(eR/np.array([MD(tt,rho) for tt in Te]))[list(Te).index(t)]),
                                     (x*1.02,1.0),fontsize=8)
ax2.set_xlabel(r'$v_p/v_{th,e}$ at birth'); ax2.set_ylabel('eRPA / MD-fit range ratio')
ax2.set_title('Divergence is purely low-velocity'); ax2.invert_xaxis(); ax2.grid(alpha=.3)
ax2.text(0.95,1.55,'models agree\nat the peak',fontsize=9,ha='left',color='0.3')
fig.tight_layout(); fig.savefig('alpha_range_models.png',dpi=140); fig.savefig('alpha_range_models.pdf')
print('wrote alpha_range_models.png/.pdf')