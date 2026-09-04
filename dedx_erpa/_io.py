# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""Self-contained readers for dedx-erpa output files (NumPy only, no FAC).

These parse the plain-text ``dedx.dat`` (electronic) and ``dedx_nuc.dat``
(electronic + nuclear/ionic total) result files written by the dedx-erpa
pipeline into NumPy arrays plus a header dictionary.  They deliberately avoid
importing the compute modules (``rd``, ``dedx``), which require FAC/``pfac``, so
the public units-aware API in :mod:`dedx_erpa` can read results in any
environment -- including continuous integration without FAC installed.
"""
from __future__ import annotations

import os
import re

import numpy as np


def _read_header(path):
    """Parse the ``# key = value`` comment header into a dict.

    Values are returned as float when possible, a list of floats for
    whitespace-separated numbers (e.g. compound ``zt``), otherwise the raw
    string (e.g. ``npot = gk``).
    """
    meta = {}
    with open(path) as fh:
        for line in fh:
            if not line.lstrip().startswith("#"):
                break
            body = line.lstrip()[1:].strip()
            if "=" not in body:
                continue
            key, _, val = body.partition("=")
            key, val = key.strip(), val.strip()
            try:
                meta[key] = float(val)
            except ValueError:
                try:
                    meta[key] = [float(p) for p in val.split()]
                except ValueError:
                    meta[key] = val
    return meta


def _column_header(path):
    """Return the tokens of the last comment line (the column header)."""
    last = None
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                last = line.lstrip()[1:].strip()
            else:
                break
    return last.split() if last else []


def read_electronic(od):
    """Read ``<od>/dedx.dat`` -> ``(meta, columns)``.

    Columns (NumPy arrays): ``energy`` [MeV/amu], ``electronic`` stopping cross
    section [1e-15 eV cm^2 per target atom], ``range`` [mg/cm^2].
    """
    path = os.path.join(od, "dedx.dat")
    meta = _read_header(path)
    data = np.loadtxt(path, comments="#", unpack=True)
    return meta, {"energy": data[0], "electronic": data[1], "range": data[2]}


def read_total(od):
    """Read ``<od>/dedx_nuc.dat`` -> ``(meta, columns)``.

    Always provides ``energy``, ``electronic``, ``total`` and ``range``; each
    nuclear/ionic potential present is keyed as ``nuclear_<name>`` (e.g.
    ``nuclear_gk``).
    """
    path = os.path.join(od, "dedx_nuc.dat")
    meta = _read_header(path)
    cols = _column_header(path)
    data = np.loadtxt(path, comments="#", unpack=True)
    out = {}
    for i, name in enumerate(cols):
        if i == 0:
            out["energy"] = data[i]
        elif name == "dEdx_e":
            out["electronic"] = data[i]
        elif name == "dEdx_tot":
            out["total"] = data[i]
        elif name == "range":
            out["range"] = data[i]
        else:
            m = re.match(r"dEdx_n\[(.+)\]", name)
            if m:
                out[f"nuclear_{m.group(1)}"] = data[i]
    return meta, out
