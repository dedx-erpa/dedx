# SPDX-License-Identifier: GPL-3.0-or-later
import os

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def coldal():
    """Path to a finished cold-Al result directory.

    Uses the committed trimmed fixture under ``tests/data/sample`` (so CI has it
    without FAC); falls back to a full local ``ColdAl/`` run if present.
    """
    fixture = os.path.join(_HERE, "data", "sample")
    full = os.path.join(os.path.dirname(_HERE), "ColdAl")
    for d in (full, fixture):
        if os.path.exists(os.path.join(d, "dedx.dat")):
            return d
    pytest.skip("no sample result directory available")
