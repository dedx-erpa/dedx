# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""dedx-erpa -- enhanced RPA/LDA ion stopping power, with PlasmaPy-friendly
(``astropy.units``) result objects.

The heavy compute pipeline (FAC average-atom densities + the Fortran ``dedx``
kernel) is driven by the command-line tool ``dedx.py``.  This package wraps the
results in :mod:`astropy.units` quantities so they interoperate cleanly with
PlasmaPy and the wider scientific-Python ecosystem.

Examples
--------
>>> import dedx_erpa
>>> res = dedx_erpa.read_total("ColdAl")     # a finished result directory
>>> res.energy.unit, res.total.unit
(Unit("MeV / u"), Unit("eV cm2"))
>>> res.total[0] > res.electronic[0]         # nuclear term adds to electronic
np.True_

Notes
-----
``energy`` is energy *per atomic mass unit* (per nucleon); multiply by the
projectile mass number A to get total kinetic energy.  ``electronic``/``total``
are stopping *cross sections* per target atom (energy x area); multiply by the
target number density to obtain dE/dx.  ``range`` is the **physical** (per-ion)
CSDA range: the pipeline applies the projectile mass-number factor internally,
so no per-nucleon correction is needed here.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import astropy.units as u

from . import _io

__version__ = "0.1.0"
__all__ = [
    "StoppingResult",
    "TotalStoppingResult",
    "__version__",
    "read_stopping",
    "read_total",
]

# Unit conventions of the dedx-erpa output files.
_E = u.MeV / u.u                  # energy per atomic mass unit (per nucleon)
_S = 1e-15 * u.eV * u.cm ** 2     # stopping cross section per target atom
_R = u.mg / u.cm ** 2             # CSDA areal-density range


@dataclass
class StoppingResult:
    """Electronic stopping result, as astropy-units quantities."""

    energy: u.Quantity            # [MeV / u]
    electronic: u.Quantity        # [eV cm^2] per target atom
    range: u.Quantity             # [mg / cm^2]
    z_target: int
    z_projectile: int
    temperature: u.Quantity       # electron temperature [eV]
    mass_density: u.Quantity      # [g / cm^3]
    mean_ionization: float


@dataclass
class TotalStoppingResult:
    """Total (electronic + nuclear/ionic) stopping result, units-aware."""

    energy: u.Quantity
    electronic: u.Quantity
    total: u.Quantity
    range: u.Quantity
    nuclear: dict = field(default_factory=dict)   # {potential_name: Quantity}
    z_projectile: int = 0
    temperature: u.Quantity = None                # electron temperature [eV]
    ion_temperature: u.Quantity = None            # ion temperature [eV]


def _as_int(value):
    return int(value[0] if isinstance(value, list) else value)


def read_stopping(od):
    """Load a finished electronic-stopping directory as a `StoppingResult`."""
    meta, c = _io.read_electronic(od)
    return StoppingResult(
        energy=c["energy"] * _E,
        electronic=c["electronic"] * _S,
        range=c["range"] * _R,
        z_target=_as_int(meta.get("zt", 0)),
        z_projectile=_as_int(meta.get("zp", 0)),
        temperature=meta.get("Te", 0.0) * u.eV,
        mass_density=meta.get("rho", 0.0) * u.g / u.cm ** 3,
        mean_ionization=float(meta.get("zbar", 0.0)),
    )


def read_total(od):
    """Load a finished total-stopping directory as a `TotalStoppingResult`."""
    meta, c = _io.read_total(od)
    nuclear = {
        k.replace("nuclear_", ""): v * _S
        for k, v in c.items()
        if k.startswith("nuclear_")
    }
    te = meta.get("Te", None)
    ti = meta.get("Ti", None)
    return TotalStoppingResult(
        energy=c["energy"] * _E,
        electronic=c["electronic"] * _S,
        total=c["total"] * _S,
        range=c["range"] * _R,
        nuclear=nuclear,
        z_projectile=_as_int(meta.get("zp", 0)),
        temperature=(te * u.eV) if te is not None else None,
        ion_temperature=(ti * u.eV) if ti is not None else None,
    )
