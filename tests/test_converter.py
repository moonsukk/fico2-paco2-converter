"""Tests for the FiCO2 <-> PaCO2 converter.

Run with:  python -m pytest        (or)   python tests/test_converter.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fico2_paco2_converter import (  # noqa: E402
    Params, params_from_baseline, params_resting_default,
    paco2_to_fico2, fico2_to_paco2, fico2_to_paco2_numeric,
    mmhg_to_kpa, kpa_to_mmhg, VA_REST, DEFAULT_BASELINE,
)


def test_note_worked_example():
    """review_derivation Example A: baseline 30, target 60 -> FiCO2 ~8.1%."""
    p = params_from_baseline(30.0)
    assert math.isclose(p.VA_base, 0.863 * 200 / 30, rel_tol=1e-9)
    fi = paco2_to_fico2(60.0, p)
    assert abs(fi - 8.135) < 0.01


def test_roundtrip_baseline40():
    """Forward then reverse must return the starting value across the range."""
    p = params_from_baseline(40.0)
    for fi in (0.0, 1.0, 2.5, 5.0, 7.0, 8.0):
        pa = fico2_to_paco2(fi, p)
        back = paco2_to_fico2(pa, p)
        assert abs(back - fi) < 1e-9


def test_closed_form_matches_numeric():
    """Exact quadratic root must match the independent bisection solver."""
    p = params_from_baseline(40.0)
    for fi in (0.0, 1.0, 3.0, 5.0, 8.0):
        assert abs(fico2_to_paco2(fi, p) - fico2_to_paco2_numeric(fi, p)) < 1e-6


def test_technical_note_table1():
    """technical_note Table 1 reference values (baseline 40, S=2.69)."""
    p = params_from_baseline(40.0)
    expected = {0.0: 40.0, 2.0: 40.8, 5.0: 45.2, 7.0: 54.0, 8.0: 60.0}
    for fi, pa_exp in expected.items():
        assert abs(fico2_to_paco2(fi, p) - pa_exp) < 0.1


def test_baseline_self_consistency():
    """With a known baseline, FiCO2 = 0 returns the baseline exactly."""
    for base in (30.0, 40.0, 45.0):
        p = params_from_baseline(base)
        assert abs(fico2_to_paco2(0.0, p) - base) < 1e-9


def test_resting_fallback_uses_4_2():
    """No-baseline fallback: VA_base = 4.2 L/min, PaCO2_base = 40 mmHg."""
    p = params_resting_default()
    assert p.VA_base == VA_REST == 4.2
    assert p.PaCO2_base == DEFAULT_BASELINE == 40.0
    # steep HCVR absorbs the small VA_base offset -> FiCO2 0 ~ 40.04 mmHg
    assert abs(fico2_to_paco2(0.0, p) - 40.0) < 0.1


def test_fallback_vs_selfconsistent_are_close():
    """5% CO2 gives ~45.2 mmHg either way (fallback vs self-consistent)."""
    pa_fb = fico2_to_paco2(5.0, params_resting_default())
    pa_sc = fico2_to_paco2(5.0, params_from_baseline(40.0))
    assert abs(pa_fb - pa_sc) < 0.1
    assert abs(pa_fb - 45.2) < 0.2


def test_unit_helpers():
    assert abs(mmhg_to_kpa(760.0) - 101.325) < 0.05
    assert abs(kpa_to_mmhg(1.0) - 7.500617) < 1e-6
    assert abs(mmhg_to_kpa(kpa_to_mmhg(5.0)) - 5.0) < 1e-5


def test_negative_fico2_raises():
    try:
        fico2_to_paco2(-1.0, params_from_baseline(40.0))
    except ValueError:
        return
    raise AssertionError("negative FiCO2 should raise ValueError")




# ---- Model B (nonlinear, duration-aware) --------------------------------
from fico2_paco2_converter import (  # noqa: E402
    ParamsB, fico2_to_paco2_B, paco2_to_fico2_B,
)


def test_modelB_reduces_to_A_in_working_range():
    """At steady state and PaCO2 <= P2, Model B == Model A."""
    pA = params_from_baseline(40.0)
    pB = ParamsB()
    for fi in (0.0, 2.0, 5.0, 7.0):
        assert abs(fico2_to_paco2_B(fi, pB) - fico2_to_paco2(fi, pA)) < 1e-6


def test_modelB_duration_raises_paco2():
    """Shorter duration -> less ventilation -> higher PaCO2."""
    pB = ParamsB()
    pa_10s = fico2_to_paco2_B(7.0, pB, t_s=10)
    pa_ss = fico2_to_paco2_B(7.0, pB, t_s=None)
    assert pa_10s > pa_ss > 40.0
    # monotone in time
    seq = [fico2_to_paco2_B(7.0, pB, t_s=t) for t in (5, 30, 120, 600)]
    assert all(a >= b for a, b in zip(seq, seq[1:]))


def test_modelB_roundtrip():
    """Reverse (iterative) then forward recovers the input at any duration."""
    pB = ParamsB()
    for fi in (5.0, 12.0, 25.0):
        for t in (30, None):
            pa = fico2_to_paco2_B(fi, pB, t_s=t)
            assert abs(paco2_to_fico2_B(pa, pB, t_s=t) - fi) < 1e-4


def test_modelB_phi_bounds():
    """Duration factor is 0<phi<=1 and ->1 as t->inf."""
    pB = ParamsB()
    assert 0 < pB.phi(1) < pB.phi(100) < 1
    assert abs(pB.phi(float('inf')) - 1.0) < 1e-12

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS  {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")
