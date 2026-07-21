#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fico2_paco2_converter.py
========================
Reversible conversion between inspired CO2 fraction (FiCO2) and arterial /
end-tidal CO2 partial pressure (PaCO2 = PetCO2), using the *simplified*
hypercapnic-ventilatory-response (HCVR) method with a single fixed slope
S = 2.69 L.min^-1.mmHg^-1  (Hirshman et al. 1975: mean 2.69 +/- 0.19 SEM, n = 44; individual range 1.00-5.95).

Assumption: PaCO2 = PetCO2 = P_ACO2 (ideal alveolar / arterial equilibration).

--------------------------------------------------------------------------
Governing equation (implicit, forward = solve for FiCO2 given PaCO2):

    PaCO2 = FiCO2*(Patm - PH2O)  +  K * VCO2 / VA(PaCO2)                (Eq. 9)

with the CO2-sensitive alveolar ventilation

    VA(PaCO2) = VA_base + S * (PaCO2 - PaCO2_base)                     (Eq. 8)

Because VA is *linear* in PaCO2, Eq. 9 is a QUADRATIC in PaCO2 and has an
exact closed-form solution -- no iterative root finding is required.

    S*PaCO2^2 + (D - S*PICO2)*PaCO2 - (PICO2*D + K*VCO2) = 0
        where  D = VA_base - S*PaCO2_base,   PICO2 = FiCO2*(Patm - PH2O)

    PaCO2 = [ (S*PICO2 - D) + sqrt((S*PICO2 - D)^2 + 4*S*(PICO2*D + K*VCO2)) ]
            / (2*S)                                                    (physical root)

Units convention for K:
    K = 0.863 when VCO2 is in mL.min^-1 (STPD) and VA is in L.min^-1 (BTPS).
    (K = (310/273)*760 / 1000 = 0.863;  see documentation.)

--------------------------------------------------------------------------
How VA_base is chosen (matches the two "currencies" seen in the literature):

  * If a baseline PaCO2 (PaCO2_base) is known, VA_base is set for
    self-consistency,  VA_base = K*VCO2/PaCO2_base, so that FiCO2 = 0 returns
    the baseline exactly.
  * If only FiCO2 is known (no baseline), we fall back to PaCO2_base = 40 mmHg
    together with a resting alveolar ventilation VA_base = 4.3 L/min -- the
    value implied by the same VCO2 = 200 and PaCO2_base = 40 (K*VCO2/40 =
    4.315 ~ 4.3), well within the normal resting range (West; Nunn's).  This is
    essentially self-consistent, so FiCO2 = 0 returns ~40.0 mmHg -- a negligible
    offset from the assumed 40.

Author: prepared for M. Kim, Isometabolism Review / Technical Note (2026).
"""

from __future__ import annotations
import argparse
import math
from dataclasses import dataclass

MMHG_PER_KPA = 7.500617       # 1 kPa = 7.500617 mmHg
KPA_PER_MMHG = 0.1333224      # 1 mmHg = 0.1333224 kPa

# Resting alveolar ventilation fallback, L/min (BTPS).  Uses 4.3 -- the value
# implied by the same VCO2 = 200 and PaCO2_base = 40 (K*VCO2/40 = 4.315 ~ 4.3),
# which sits within the normal resting range (West; Nunn's ~4-5 L/min) and is
# essentially self-consistent, so FiCO2 = 0 returns ~40.0 mmHg.
VA_REST = 4.3

DEFAULT_BASELINE = 40.0        # normocapnic resting PaCO2, mmHg


@dataclass
class Params:
    """Physiological / environmental constants for the simplified converter."""
    S: float = 2.69            # HCVR slope, L.min^-1.mmHg^-1 (Hirshman 1975)
    VCO2: float = 200.0        # CO2 output, mL.min^-1 STPD (Nunn's: 150-200)
    PaCO2_base: float = 40.0   # baseline (normocapnic) PaCO2, mmHg
    VA_base: float | None = None  # baseline alveolar ventilation, L.min^-1 BTPS
    Patm: float = 760.0        # barometric pressure, mmHg
    PH2O: float = 47.0         # saturated water-vapour pressure at 37 C, mmHg
    K: float = 0.863           # STPD->BTPS / fraction->pressure constant

    def __post_init__(self):
        # If VA_base is not supplied, fix it so that FiCO2 = 0 reproduces the
        # baseline PaCO2 exactly (self-consistency):  VA_base = K*VCO2/PaCO2_base
        if self.VA_base is None:
            self.VA_base = self.K * self.VCO2 / self.PaCO2_base

    @property
    def Pdry(self) -> float:
        """Dry-gas partial-pressure scaling (Patm - PH2O), mmHg."""
        return self.Patm - self.PH2O

    @property
    def D(self) -> float:
        """Constant part of the ventilation denominator: VA_base - S*PaCO2_base."""
        return self.VA_base - self.S * self.PaCO2_base

    def VA(self, paco2: float) -> float:
        """CO2-sensitive alveolar ventilation at a given PaCO2 (L.min^-1)."""
        return self.VA_base + self.S * (paco2 - self.PaCO2_base)


# --------------------------------------------------------------------------
# Parameter builders for the two "known-inputs" scenarios
# --------------------------------------------------------------------------
def params_from_baseline(paco2_base: float, **kw) -> Params:
    """Baseline PaCO2 is known -> VA_base is set self-consistently.

    VA_base = K*VCO2/PaCO2_base, so FiCO2 = 0 reproduces PaCO2_base exactly.
    Use for both directions whenever the user supplies a baseline PaCO2.
    """
    return Params(PaCO2_base=paco2_base, VA_base=None, **kw)


def params_resting_default(**kw) -> Params:
    """Only FiCO2 is known (no baseline) -> literature resting fallback.

    PaCO2_base = 40 mmHg (HCVR reference point) and VA_base = 4.3 L/min
    (the self-consistent value K*VCO2/40 ~ 4.3, within the normal resting
    range). These are essentially self-consistent, so FiCO2 = 0 returns ~40.0 mmHg.
    """
    return Params(PaCO2_base=DEFAULT_BASELINE, VA_base=VA_REST, **kw)


# --------------------------------------------------------------------------
# Forward direction:  PaCO2  ->  FiCO2   (explicit, algebraic)
# --------------------------------------------------------------------------
def paco2_to_fico2(paco2_mmhg: float, p: Params | None = None,
                   as_percent: bool = True) -> float:
    """Given a target PaCO2 (mmHg), return the inspired CO2 fraction FiCO2.

    FiCO2 = [ PaCO2 - K*VCO2 / VA(PaCO2) ] / (Patm - PH2O)
    Returns a percentage by default, or a fraction if as_percent=False.
    """
    p = p or Params()
    va = p.VA(paco2_mmhg)
    if va <= 0:
        raise ValueError(f"Non-physical: alveolar ventilation <= 0 at PaCO2={paco2_mmhg}")
    pico2 = paco2_mmhg - p.K * p.VCO2 / va      # inspired CO2 partial pressure, mmHg
    fico2 = pico2 / p.Pdry
    return fico2 * 100.0 if as_percent else fico2


# --------------------------------------------------------------------------
# Reverse direction:  FiCO2  ->  PaCO2   (exact closed-form quadratic root)
# --------------------------------------------------------------------------
def fico2_to_paco2(fico2, p: Params | None = None,
                   is_percent: bool = True) -> float:
    """Given inspired CO2 (percent by default), return PaCO2 in mmHg.

    Solves the quadratic exactly and returns the physical (upper) root.
    """
    p = p or Params()
    fico2_frac = fico2 / 100.0 if is_percent else float(fico2)
    if fico2_frac < 0:
        raise ValueError("FiCO2 cannot be negative.")
    pico2 = fico2_frac * p.Pdry                  # inspired CO2 partial pressure, mmHg
    a = p.S
    b = p.D - p.S * pico2
    c = -(pico2 * p.D + p.K * p.VCO2)
    disc = b * b - 4 * a * c
    if disc < 0:
        raise ValueError("No real solution (negative discriminant).")
    paco2 = (-b + math.sqrt(disc)) / (2 * a)     # '+' root is the physical one
    if p.VA(paco2) <= 0:
        raise ValueError("Solved PaCO2 lies below the apnoeic threshold (VA<=0).")
    return paco2


# --------------------------------------------------------------------------
# Independent numerical cross-check (root finder) -- for verification only
# --------------------------------------------------------------------------
def fico2_to_paco2_numeric(fico2, p: Params | None = None,
                           is_percent: bool = True) -> float:
    """Bisection solver bracketed to the physical branch (VA>0). Used to
    confirm the closed-form result; not needed for normal use."""
    p = p or Params()
    fico2_frac = fico2 / 100.0 if is_percent else float(fico2)
    pico2 = fico2_frac * p.Pdry

    def resid(x):
        return x - (pico2 + p.K * p.VCO2 / p.VA(x))

    lo = p.PaCO2_base - p.VA_base / p.S + 1e-6   # apnoeic threshold (VA -> 0)
    hi = 400.0
    flo, fhi = resid(lo), resid(hi)
    if flo * fhi > 0:
        raise ValueError("Root not bracketed on the physical branch.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = resid(mid)
        if abs(fm) < 1e-12:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# ==========================================================================
# TIME-DEPENDENT MODEL  (duration as an explicit variable)
# ==========================================================================
# The steady-state model (above) assumes the ventilatory response to the
# imposed CO2 has fully developed.  The time-dependent model adds ONE extra
# variable -- the duration t (seconds) of the hypercapnic challenge -- because
# the response is not instantaneous.
#
# The acute response develops through a fast (peripheral, carotid-body) and a
# slower central (medullary) component.  FIVE studies report both time
# constants in normoxic humans during step hypercapnia:
#
#   study                        n    tau_fast     tau_central    λ (peripheral share)
#   ---------------------------  ---  -----------  -------------  --------------------
#   Swanson & Bellville (1975)    1   17.5 s       75.0 s         ~0.50
#   Ward & Bellville   (1983)     6    5.7 s        72.0 s         0.29   (g2/(g1+g2))
#   Bellville et al.   (1979)     7   14.8+-11.1s  180.1+-96.0 s  0.34   (g2/(g1+g2))
#   Dahan et al.       (1990)     9    9.8+- 3.5s  146.6+-48.8 s  0.30   (Gp/(Gc+Gp))
#   Berkenbosch et al. (1992)    10    7.6 s       127.5 s        0.26   (Sp/(Sc+Sp))
#
# Weighting each study by its number of subjects gives the defaults used here:
#
#       tau_fast    = 9.7 s       tau_central = 132.2 s      f = 0.295
#
# Both lie inside the ranges of 8-26 s and 65-180 s summarised by Tansley et
# al. (1998), and tau_fast agrees with the ~15 s adopted by Cunningham et al.
# (1986).  The five source studies used comparable CO2 doses: their end-tidal
# steps of +7 to +11 mmHg correspond, via the steady-state model, to an
# equivalent inspired fraction of ~5-7% CO2.
#
# The fraction of the acute steady-state response present at time t is
#
#       phi(t) = λ (1 - e^{-t/tau_fast}) + (1 - λ)(1 - e^{-t/tau_central}),
#
# rising from 0 to 1 over the first ~5-10 min.  (Over hours-days the response
# acclimatises, but that regime is variable between studies and is NOT
# modelled here.)
#
# The EFFECTIVE slope at time t is  S_eff(t) = phi(t) * S, and ventilation is
#   VA(PaCO2, t) = VA_base + phi(t) * S * (PaCO2 - PaCO2_base).
# For a FIXED t this is linear in PaCO2, so the time-dependent model is simply
# the steady-state model evaluated with the reduced slope S_eff(t): early in a
# challenge S_eff is small, less CO2 is offloaded, and PaCO2 is transiently
# higher.  It converges to the steady-state model as t -> inf.
#
# The steady-state slope S itself varies widely between individuals and method
# (reported ~1-6 L/min/mmHg; Read 1967; Hirshman 1975 give 1.00-5.95); the
# default 2.69 is the population mean and is adjustable.  The model is intended
# for the human working range: PaCO2 <= cap_paco2 (~80 mmHg); above that, use a
# measured PaCO2.


@dataclass
class ParamsB:
    """Parameters for the time-dependent converter."""
    S: float = 2.69            # steady-state HCVR slope, L/min/mmHg (reported 1.00-5.95)
    VCO2: float = 200.0        # CO2 output, mL/min STPD
    PaCO2_base: float = 40.0   # normocapnic baseline, mmHg
    VA_base: float | None = None
    Patm: float = 760.0
    PH2O: float = 47.0
    K: float = 0.863
    tau_fast: float = 9.7      # peripheral time constant, s  (subject-weighted mean of Swanson
                               #   1975 17.5 / Ward-Bellville 1983 5.7 / Bellville 1979 14.8 /
                               #   Dahan 1990 9.8 / Berkenbosch 1992 7.6; n = 33)
    tau_central: float = 132.2  # central time constant, s   (subject-weighted mean of Swanson
                               #   1975 75 / Ward-Bellville 1983 72 / Bellville 1979 180.1 /
                               #   Dahan 1990 146.6 / Berkenbosch 1992 127.5; n = 33)
    frac_fast: float = 0.295   # peripheral share of the steady-state response (subject-weighted
                               #   mean of Bellville 1979 0.34, Dahan 1990 0.30, Ward-Bellville
                               #   1983 0.29, Berkenbosch 1992 0.26; n = 32)
    cap_paco2: float = 80.0    # upper validity bound for the human model, mmHg

    def __post_init__(self):
        if self.VA_base is None:
            self.VA_base = self.K * self.VCO2 / self.PaCO2_base

    @property
    def Pdry(self) -> float:
        return self.Patm - self.PH2O

    def phi(self, t_s: float | None) -> float:
        """Fraction of the acute steady-state response present at time t_s
        seconds after onset.  t_s = None or inf -> 1 (steady state)."""
        if t_s is None or math.isinf(t_s):
            return 1.0
        return (self.frac_fast * (1.0 - math.exp(-t_s / self.tau_fast))
                + (1.0 - self.frac_fast) * (1.0 - math.exp(-t_s / self.tau_central)))

    def S_eff(self, t_s: float | None) -> float:
        """Effective HCVR slope at time t_s:  S_eff(t) = phi(t) * S."""
        return self.phi(t_s) * self.S

    def VA(self, paco2: float, t_s: float | None) -> float:
        """Alveolar ventilation present at PaCO2 and time t_s."""
        return self.VA_base + self.S_eff(t_s) * (paco2 - self.PaCO2_base)


def paco2_to_fico2_B(paco2_mmhg: float, p: ParamsB | None = None,
                     t_s: float | None = None, as_percent: bool = True) -> float:
    """Time-dependent forward (explicit): PaCO2 + duration t_s (s) -> FiCO2."""
    p = p or ParamsB()
    va = p.VA(paco2_mmhg, t_s)
    if va <= 0:
        raise ValueError(f"Non-physical: VA <= 0 at PaCO2={paco2_mmhg}")
    pico2 = paco2_mmhg - p.K * p.VCO2 / va
    fico2 = pico2 / p.Pdry
    return fico2 * 100.0 if as_percent else fico2


def fico2_to_paco2_B(fico2, p: ParamsB | None = None, t_s: float | None = None,
                     is_percent: bool = True) -> float:
    """Time-dependent reverse: FiCO2 + duration t_s (seconds) -> PaCO2.

    For a fixed t the model is linear in PaCO2 (slope S_eff(t) = phi(t)*S), so
    this is solved by bracketed bisection on the monotone residual -- robust
    even in the t -> 0 limit where S_eff -> 0 (no ventilatory response yet).
    """
    p = p or ParamsB()
    frac = fico2 / 100.0 if is_percent else float(fico2)
    if frac < 0:
        raise ValueError("FiCO2 cannot be negative.")
    pico2 = frac * p.Pdry

    def resid(paco2: float) -> float:
        return paco2 - pico2 - p.K * p.VCO2 / p.VA(paco2, t_s)

    lo, hi = p.PaCO2_base, 600.0
    if resid(lo) > 0:
        return lo
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        fm = resid(mid)
        if abs(fm) < 1e-9:
            return mid
        if fm < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------
# Convenience unit helpers
# --------------------------------------------------------------------------
def mmhg_to_kpa(x: float) -> float:
    return x * KPA_PER_MMHG


def kpa_to_mmhg(x: float) -> float:
    return x * MMHG_PER_KPA


# --------------------------------------------------------------------------
# Interactive, input-driven front end
# --------------------------------------------------------------------------
def _ask_float(prompt: str, default: float | None = None) -> float:
    """Prompt until the user types a valid number. Blank -> default (if given)."""
    while True:
        suffix = f" [{default}]" if default is not None else ""
        raw = input(f"{prompt}{suffix}: ").strip()
        if not raw and default is not None:
            return default
        try:
            return float(raw)
        except ValueError:
            print("  Please enter a number (e.g. 40 or 5.5).")


def _ask_yes_no(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no. Blank -> default."""
    tag = "Y/n" if default else "y/N"
    while True:
        raw = input(f"{prompt} [{tag}]: ").strip().lower()
        if not raw:
            return default
        if raw in ("y", "yes"):
            return True
        if raw in ("n", "no"):
            return False
        print("  Please answer y or n.")


def _cap_note(paco2: float) -> str:
    return ("   [!] PaCO2 > 80 mmHg -- outside the human working range; "
            "use a measured PaCO2." if paco2 > 80.0 else "")


def run_paco2_to_fico2(t_s: float | None = None) -> None:
    """Direction 1: target PaCO2 (+ baseline) -> FiCO2.  t_s given => time-dependent."""
    print("\n-- Estimate FiCO2 from a target PaCO2 --")
    paco2 = _ask_float("Target PaCO2 (the hypercapnic level), mmHg")
    paco2_base = _ask_float("Baseline (normocapnic) PaCO2, mmHg", default=DEFAULT_BASELINE)
    if t_s is None:
        p = params_from_baseline(paco2_base)
        fi, va, seff = paco2_to_fico2(paco2, p), p.VA(paco2), p.S
    else:
        p = ParamsB(PaCO2_base=paco2_base)
        fi, va, seff = paco2_to_fico2_B(paco2, p, t_s=t_s), p.VA(paco2, t_s), p.S_eff(t_s)
    print(f"\n  PaCO2 = {paco2:.2f} mmHg ({mmhg_to_kpa(paco2):.3f} kPa), baseline {paco2_base:.0f}")
    print(f"  effective slope S_eff = {seff:.2f} L/min/mmHg ;  VA = {va:.2f} L/min")
    print(f"  ->  FiCO2 = {fi:.3f} %  (PICO2 = {fi/100*(p.Patm-p.PH2O):.1f} mmHg)")
    if _cap_note(paco2):
        print(_cap_note(paco2))


def run_fico2_to_paco2(t_s: float | None = None) -> None:
    """Direction 2: FiCO2 (+ baseline) -> PaCO2.  t_s given => time-dependent."""
    print("\n-- Estimate PaCO2 from a given FiCO2 --")
    fico2 = _ask_float("Inspired CO2, FiCO2 (percent, e.g. 5)")
    knows_base = _ask_yes_no("Do you know the subject's baseline (normocapnic) PaCO2?",
                             default=False)
    base = _ask_float("Baseline PaCO2, mmHg", default=DEFAULT_BASELINE) if knows_base else DEFAULT_BASELINE
    if t_s is None:
        p = params_from_baseline(base) if knows_base else params_resting_default()
        pa, seff = fico2_to_paco2(fico2, p), p.S
        va = p.VA(pa)
    else:
        p = ParamsB(PaCO2_base=base)
        if not knows_base:
            p.VA_base = VA_REST      # literature resting VA when baseline unknown
        pa, seff = fico2_to_paco2_B(fico2, p, t_s=t_s), p.S_eff(t_s)
        va = p.VA(pa, t_s)
    print(f"\n  FiCO2 = {fico2:.3f} %  (PICO2 = {fico2/100*(p.Patm-p.PH2O):.1f} mmHg)")
    print(f"  effective slope S_eff = {seff:.2f} L/min/mmHg ;  VA = {va:.2f} L/min")
    print(f"  ->  PaCO2 = {pa:.2f} mmHg = {mmhg_to_kpa(pa):.3f} kPa  "
          f"(rise of {pa - base:+.2f} mmHg from baseline)")
    if _cap_note(pa):
        print(_cap_note(pa))


def interactive() -> None:
    """Top-level interactive menu (choose the model, then a direction)."""
    print("=" * 68)
    print(" FiCO2 <-> PaCO2 converter")
    print("=" * 68)
    print(" Which model?")
    print("   S) steady-state    -- response fully developed; closed form (default S=2.69)")
    print("   T) time-dependent  -- you also give the challenge DURATION; uses S_eff(t)=phi(t)*S")
    t_s = None
    while True:
        m = input(" Enter S or T (q to quit): ").strip().lower()
        if m in ("q", "quit", "exit"):
            return
        if m in ("s", "t", "a", "b"):   # a/b accepted for backwards compatibility
            break
        print("  Please type S or T.")
    if m in ("t", "b"):
        mins = _ask_float("Challenge duration so far, minutes", default=5.0)
        t_s = mins * 60.0
        pb = ParamsB()
        print(f"# time-dependent: at t = {mins:g} min, phi = {pb.phi(t_s):.3f} -> "
              f"S_eff = {pb.S_eff(t_s):.2f} L/min/mmHg (steady-state S = {pb.S})")

    print("\n What do you know / want?")
    print("   1) I have a target PaCO2  ->  estimate FiCO2")
    print("   2) I have an inspired FiCO2  ->  estimate PaCO2")
    while True:
        choice = input(" Enter 1 or 2 (q to quit): ").strip().lower()
        if choice == "1":
            run_paco2_to_fico2(t_s)
            return
        if choice == "2":
            run_fico2_to_paco2(t_s)
            return
        if choice in ("q", "quit", "exit"):
            return
        print("  Please type 1, 2, or q.")


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------
def _selftest() -> None:
    """Reproduce the worked examples and check both directions."""
    print("=" * 64)
    print("SELF-TEST")
    print("=" * 64)

    # (A) Note's worked example: baseline PaCO2 = 30 -> target 60 -> FiCO2 ~8.1%.
    note = params_from_baseline(30.0)  # VA_base auto = 5.753 L/min
    f = paco2_to_fico2(60.0, note)
    p = fico2_to_paco2(f, note)
    print(f"[note] PaCO2 60 -> FiCO2 {f:6.3f}%  (expected ~8.1%)")
    print(f"[note] FiCO2 {f:.3f}% -> PaCO2 {p:7.4f}  (expected 60)")
    assert abs(f - 8.135) < 0.01
    assert abs(p - 60.0) < 1e-6

    # (B) Default normocapnic baseline (PaCO2 = 40), round-trip + numeric check.
    dp = params_from_baseline(40.0)  # VA_base auto = 4.315 L/min
    print(f"\n[baseline=40] VA_base = {dp.VA_base:.3f} L/min (self-consistent)")
    print(f"{'FiCO2 %':>8} {'PaCO2 mmHg':>12} {'PaCO2 kPa':>11} {'numeric':>10} {'round-trip %':>13}")
    for fi in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0):
        pa = fico2_to_paco2(fi, dp)
        num = fico2_to_paco2_numeric(fi, dp)
        rt = paco2_to_fico2(pa, dp)
        assert abs(pa - num) < 1e-6, "closed-form vs numeric mismatch"
        assert abs(rt - fi) < 1e-9, "round-trip mismatch"
        print(f"{fi:8.1f} {pa:12.3f} {mmhg_to_kpa(pa):11.3f} {num:10.3f} {rt:13.4f}")

    # (C) Resting fallback (no baseline): VA_base = 4.2, PaCO2_base = 40.
    rp = params_resting_default()
    pa0 = fico2_to_paco2(0.0, rp)
    print(f"\n[resting fallback] VA_base = {rp.VA_base:.1f} L/min, "
          f"FiCO2 0% -> PaCO2 {pa0:.3f} mmHg (expect ~40.04, <0.1 mmHg offset)")
    assert abs(pa0 - 40.0) < 0.1
    print("\nAll self-tests passed.")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Reversible FiCO2 <-> PaCO2 converter (simplified HCVR, S=2.69).")
    ap.add_argument("--fico2", type=float, help="inspired CO2 (percent) -> PaCO2")
    ap.add_argument("--paco2", type=float, help="target PaCO2 (mmHg) -> FiCO2")
    ap.add_argument("--paco2-kpa", type=float, help="target PaCO2 (kPa) -> FiCO2")
    ap.add_argument("--baseline", type=float, default=None,
                    help="baseline PaCO2, mmHg. If omitted with --fico2, uses the "
                         "resting fallback (PaCO2_base=40, VA_base=4.2 L/min).")
    ap.add_argument("--vco2", type=float, default=200.0,
                    help="CO2 output, mL/min STPD (default 200)")
    ap.add_argument("--va-base", type=float, default=None,
                    help="override baseline alveolar ventilation, L/min")
    ap.add_argument("--slope", type=float, default=2.69,
                    help="HCVR slope S, L/min/mmHg (default 2.69)")
    ap.add_argument("--patm", type=float, default=760.0, help="barometric pressure, mmHg")
    ap.add_argument("--selftest", action="store_true", help="run built-in verification")
    ap.add_argument("--interactive", "-i", action="store_true",
                    help="prompt for values interactively")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return

    # No conversion on the command line -> interactive menu.
    if args.interactive or not (args.fico2 is not None or args.paco2 is not None
                                or args.paco2_kpa is not None):
        interactive()
        return

    # --- FiCO2 -> PaCO2 (respects the baseline/resting-fallback rule) ---
    if args.fico2 is not None:
        if args.baseline is not None or args.va_base is not None:
            base = args.baseline if args.baseline is not None else DEFAULT_BASELINE
            p = Params(S=args.slope, VCO2=args.vco2, PaCO2_base=base,
                       VA_base=args.va_base, Patm=args.patm)
        else:
            p = params_resting_default(S=args.slope, VCO2=args.vco2, Patm=args.patm)
        pa = fico2_to_paco2(args.fico2, p)
        print(f"# PaCO2_base={p.PaCO2_base:.1f} mmHg, VA_base={p.VA_base:.3f} L/min, "
              f"S={p.S}, VCO2={p.VCO2}")
        print(f"FiCO2 = {args.fico2:.3f} %  ->  PaCO2 = {pa:.2f} mmHg = "
              f"{mmhg_to_kpa(pa):.3f} kPa")

    # --- PaCO2 -> FiCO2 (baseline required; defaults to 40) ---
    if args.paco2 is not None or args.paco2_kpa is not None:
        base = args.baseline if args.baseline is not None else DEFAULT_BASELINE
        p = Params(S=args.slope, VCO2=args.vco2, PaCO2_base=base,
                   VA_base=args.va_base, Patm=args.patm)
        print(f"# PaCO2_base={p.PaCO2_base:.1f} mmHg, VA_base={p.VA_base:.3f} L/min, "
              f"S={p.S}, VCO2={p.VCO2}")
        if args.paco2 is not None:
            fi = paco2_to_fico2(args.paco2, p)
            print(f"PaCO2 = {args.paco2:.2f} mmHg  ->  FiCO2 = {fi:.3f} %")
        if args.paco2_kpa is not None:
            pa_mmhg = kpa_to_mmhg(args.paco2_kpa)
            fi = paco2_to_fico2(pa_mmhg, p)
            print(f"PaCO2 = {args.paco2_kpa:.3f} kPa ({pa_mmhg:.2f} mmHg)  ->  "
                  f"FiCO2 = {fi:.3f} %")


if __name__ == "__main__":
    main()


# --------------------------------------------------------------------------
# Preferred names (the note calls these the steady-state and time-dependent
# models).  The older Model A / Model B names are kept as aliases.
# --------------------------------------------------------------------------
ParamsSteadyState = Params
ParamsTimeDependent = ParamsB
paco2_to_fico2_timedep = paco2_to_fico2_B
fico2_to_paco2_timedep = fico2_to_paco2_B
