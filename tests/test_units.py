# SPDX-License-Identifier: GPL-3.0-or-later
"""Units-API tests (require astropy; the PlasmaPy interoperability surface)."""
import numpy as np
import pytest

pytest.importorskip("astropy")
import astropy.units as u

import dedx_erpa


def test_read_stopping_units(coldal):
    r = dedx_erpa.read_stopping(coldal)
    assert r.energy.unit == u.MeV / u.u
    assert r.electronic.unit.is_equivalent(u.eV * u.cm ** 2)
    assert r.range.unit.is_equivalent(u.mg / u.cm ** 2)
    assert r.z_target == 13 and r.z_projectile == 1
    assert r.temperature.to_value(u.eV) > 0
    assert r.mass_density.to_value(u.g / u.cm ** 3) > 0
    assert r.mean_ionization > 0


def test_read_total_units_and_additivity(coldal):
    r = dedx_erpa.read_total(coldal)
    assert r.total.unit.is_equivalent(u.eV * u.cm ** 2)
    assert "gk" in r.nuclear
    nuc = sum(v for v in r.nuclear.values())
    assert np.allclose(
        r.total.to_value(u.eV * u.cm ** 2),
        (r.electronic + nuc).to_value(u.eV * u.cm ** 2),
        rtol=1e-3,
    )


def test_version():
    assert dedx_erpa.__version__
