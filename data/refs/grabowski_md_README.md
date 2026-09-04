# Grabowski (2013) classical-MD stopping — digitization slot

Ref: P. E. Grabowski, M. P. Surh, D. F. Richards, F. R. Graziani, M. S. Murillo,
"Molecular Dynamics Simulations of Classical Stopping Power," PRL 111, 215002 (2013).

To complete `validation/grabowski_compare.py`, digitize the MD stopping points
(WebPlotDigitizer, same workflow as the Frenje Fig. 5 digitization) into:

  grabowski_md_g03.csv   (Gamma_e ~ 0.3 panel)
  grabowski_md_g10.csv   (Gamma_e ~ 1  panel)
  grabowski_md_g30.csv   (Gamma_e ~ 3  panel)

Two columns, comma-separated, no header:
  col 1 = v_p / v_th,e   (projectile velocity / electron thermal velocity)
  col 2 = dE/dx in the SAME units as the plotted e-RPA curve
          [1e-15 eV cm^2 / atom].  If Grabowski is in reduced units, convert, or
          switch the y-axis of grabowski_compare.py to the reduced form and
          convert the e-RPA curve instead (one scale factor per panel).

Match the panel Gamma to Grabowski's simulated coupling, and record their exact
(n, T, mass ratio, projectile) so the e-RPA CONDITIONS list can be retuned to the
same Gamma. NOTE the classical-vs-degenerate caveat in the script docstring: the
comparison is strict only where Theta > 1 (the g03 panel); the g10/g30 panels are
indicative because the e-RPA plasma is degenerate there while Grabowski's is classical.
