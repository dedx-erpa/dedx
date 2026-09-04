# SPDX-License-Identifier: GPL-3.0-or-later
"""Parser tests (NumPy only -- run without FAC or astropy)."""
import numpy as np

from dedx_erpa import _io


def test_read_electronic(coldal):
    meta, c = _io.read_electronic(coldal)
    zt = meta["zt"][0] if isinstance(meta["zt"], list) else meta["zt"]
    assert int(zt) == 13            # aluminum
    assert int(meta["zp"]) == 1     # proton
    assert c["energy"][0] < c["energy"][-1]          # ascending energy grid
    assert np.all(c["electronic"] > 0)
    assert np.all(np.diff(c["range"]) > 0)           # range grows with energy
    assert len(c["energy"]) == len(c["electronic"]) == len(c["range"])


def test_read_total_additivity(coldal):
    meta, c = _io.read_total(coldal)
    assert {"energy", "electronic", "total", "range"} <= set(c)
    assert meta.get("npot") == "gk"
    nuc = sum(v for k, v in c.items() if k.startswith("nuclear_"))
    # total column must equal electronic + sum of nuclear potentials
    assert np.allclose(c["total"], c["electronic"] + nuc, rtol=1e-3)
