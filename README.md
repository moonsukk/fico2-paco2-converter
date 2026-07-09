# FiCO₂ ⇌ PaCO₂ Converter

A small, reversible converter between **inspired CO₂ fraction** (FiCO₂, e.g. "5% CO₂") and
**arterial / end-tidal CO₂ partial pressure** (PaCO₂ ≈ PetCO₂, in mmHg or kPa), for
harmonising the "dose" of CO₂ challenges across studies.

Built as a supplementary tool for a review on the **isometabolism assumption** (no change
in CMRO₂ during hypercapnia) in BOLD-MRI studies. CO₂ challenges are reported in two
incompatible currencies — inspired fraction vs. partial pressure — and no compact,
closed-form conversion existed because the classical alveolar equations set inspired CO₂
to zero and treat ventilation as fixed. This tool makes the ventilatory feedback explicit.

## The idea in one line

Raising inspired CO₂ raises PaCO₂, **but the resulting hypercapnia stimulates ventilation**,
which offloads part of the added CO₂. That feedback is the near-linear **hypercapnic
ventilatory response (HCVR)**, with slope `S`:

```
PaCO2 = FiCO2·(Patm − PH2O) + K·VCO2 / VA(PaCO2)          (mass balance with inspired CO2)
VA(PaCO2) = VA_base + S·(PaCO2 − PaCO2_base)              (linear HCVR)
```

Because `VA` is linear in PaCO₂, the reverse direction (FiCO₂ → PaCO₂) is a **quadratic
with an exact closed-form root** — no iterative solver needed.

## Two models

**Model A — steady-state (default).** Assumes the ventilatory response has fully developed.
Linear HCVR, reversible in closed form. Use for typical human challenges held long enough to
equilibrate (minutes).

**Model B — time-resolved.** Adds one variable, the challenge **duration `t`**, because the
response is not instantaneous — it develops with a fast (peripheral, τ ≈ 15 s) and a slow
(central, τ ≈ 130 s) component, so the *effective* slope grows over the first ~5–10 min:

```
φ(t)     = f·(1 − e^(−t/τ1)) + (1 − f)·(1 − e^(−t/τ2))   (f ≈ 0.20, τ1 ≈ 15 s, τ2 ≈ 130 s)
S_eff(t) = φ(t)·S      →   VA(PaCO2, t) = VA_base + S_eff(t)·(PaCO2 − PaCO2_base)
```

For a fixed `t` this is still linear in PaCO₂ (Model A with `S → S_eff(t)`), so it too is
closed-form (solved by a robust bisection that also covers the `t → 0` limit). Early in a
challenge PaCO₂ is transiently higher; by ~10 min **Model B → Model A**. Time constants from
Cunningham et al. (1986), Swanson & Bellville (1975), Bellville et al. (1979), Dahan et al.
(1990) and Tansley et al. (1998). Over hours–days the response acclimatises variably and is
**not** modelled.

**Validity.** Both models are for the human working range, **PaCO₂ ≤ ~80 mmHg**; above that the
HCVR saturates and no slope is measured — use a directly measured PaCO₂. The slope `S` varies
widely between individuals and methods (~1–6 L·min⁻¹·mmHg⁻¹) and is adjustable (default 2.69).

## Constants (defaults)

| Quantity | Value | Source |
|---|---|---|
| `K` | 0.863 | STPD→BTPS + fraction→pressure unit constant (West, 2012) |
| `Patm` | 760 mmHg | sea-level barometric pressure |
| `PH2O` | 47 mmHg | saturated vapour pressure at 37 °C |
| `VCO2` | 200 mL·min⁻¹ | resting CO₂ output (range 150–200) |
| `PaCO2_base` | 40 mmHg | normocapnic resting baseline |
| `VA_base` (self-consistent) | ≈ 4.3 L·min⁻¹ | `K·VCO2 / PaCO2_base` |
| `VA_REST` (fallback) | **4.2 L·min⁻¹** | literature resting alveolar ventilation (West; Nunn's) |
| `S` | 2.69 L·min⁻¹·mmHg⁻¹ | mean HCVR slope (Hirshman et al., 1975; range 1.16–5.95) |

## How the baseline / V̇A_base is chosen

- **Baseline PaCO₂ known** → `VA_base = K·VCO2 / PaCO2_base` (self-consistent: FiCO₂ = 0
  returns the baseline exactly).
- **Only FiCO₂ known** → assume `PaCO2_base = 40` mmHg and use the literature resting value
  `VA_base = 4.2 L·min⁻¹`. The steep HCVR absorbs the small mismatch, so FiCO₂ = 0 returns
  ≈ 40.04 mmHg (a negligible < 0.1 mmHg offset).

## Usage

### Interactive (prompts you for values)

```bash
python fico2_paco2_converter.py -i
```

You first pick **Model A or B** (B also asks for the challenge **duration**), then a direction:

- **Direction 1 (PaCO₂ → FiCO₂):** enter target PaCO₂ and baseline PaCO₂.
- **Direction 2 (FiCO₂ → PaCO₂):** enter FiCO₂; you're asked whether you know the baseline
  PaCO₂. If not, the resting fallback (`VA_base = 4.2`, `PaCO2_base = 40`) is used.

Model B prints the effective slope `S_eff(t)` and warns if the result exceeds ~80 mmHg.

### Command line (one-shot)

```bash
# FiCO2 -> PaCO2, no known baseline (resting fallback: VA_base = 4.2)
python fico2_paco2_converter.py --fico2 5

# FiCO2 -> PaCO2 with a known baseline (self-consistent VA_base)
python fico2_paco2_converter.py --fico2 5 --baseline 40

# PaCO2 -> FiCO2 (baseline defaults to 40 if omitted)
python fico2_paco2_converter.py --paco2 50 --baseline 40
python fico2_paco2_converter.py --paco2-kpa 6.7 --baseline 40

# override the slope or CO2 output
python fico2_paco2_converter.py --fico2 5 --baseline 40 --slope 4.0 --vco2 180
```

### As a library

```python
from fico2_paco2_converter import (
    params_from_baseline, params_resting_default,
    paco2_to_fico2, fico2_to_paco2,
)

# baseline known
p = params_from_baseline(40.0)
paco2_to_fico2(50.0, p)        # -> 6.24 %
fico2_to_paco2(5.0, p)         # -> 45.15 mmHg

# only FiCO2 known -> resting fallback (VA_base = 4.2 L/min)
fico2_to_paco2(5.0, params_resting_default())   # -> 45.17 mmHg
```

**Model B (time-resolved):**

```python
from fico2_paco2_converter import ParamsB, fico2_to_paco2_B, paco2_to_fico2_B

pB = ParamsB()                             # S=2.69, τ1=15s, τ2=130s, f=0.20, cap=80
fico2_to_paco2_B(7.0, pB, t_s=10)          # 7% for 10 s  -> 62.4 mmHg (transient)
fico2_to_paco2_B(7.0, pB, t_s=None)        # steady state -> 54.0 mmHg (= Model A)
pB.S_eff(60)                               # effective slope at 1 min -> 1.32
fico2_to_paco2_B(5.0, ParamsB(S=4.0))      # override the slope
```

### Jupyter notebook

`FiCO2_PaCO2_conversion.ipynb` is a step-by-step tutorial: the derivation, both directions,
the slope-`S` sensitivity, the **time-resolved Model B** with duration figures, and a
validation against reported challenges. Full detail is in the companion technical note.

## Reference values (baseline 40 mmHg, S = 2.69)

| FiCO₂ (%) | PICO₂ (mmHg) | PaCO₂ (mmHg) | PaCO₂ (kPa) |
|---|---|---|---|
| 0 | 0 | 40.0 | 5.33 |
| 2 | 14.3 | 40.8 | 5.44 |
| 5 | 35.7 | 45.2 | 6.02 |
| 7 | 49.9 | 54.0 | 7.20 |
| 8 | 57.0 | 60.0 | 8.00 |

## Tests

```bash
python tests/test_converter.py       # plain runner, no deps
python -m pytest                      # if pytest is installed
python fico2_paco2_converter.py --selftest
```

## Assumptions & limitations

- Ideal alveolar–arterial equilibration (PaCO₂ ≈ PACO₂ ≈ PetCO₂; healthy lungs).
- HCVR treated as linear over the working range (≈ 40–80 mmHg); **do not apply above ~80 mmHg**
  (the response saturates and no slope is measured — use a directly measured PaCO₂).
- A single default slope `S = 2.69` (adjustable; reported range ~1–6) — it will not match any
  individual, and varies with age, sex, body size, and method.
- Model B's time constants are population means from the literature; the hours–days
  acclimatisation regime is variable and is not modelled.
- Not for anaesthetised subjects or non-human species without re-parameterisation (see the
  technical note).

Best used to **estimate and harmonise reported CO₂-challenge magnitudes across studies**,
not to predict an individual subject's PaCO₂.

## Key references

- Hirshman, McCullough & Weil (1975). Normal values for hypoxic and hypercapnic ventilatory
  drives in man. *J Appl Physiol* 38(6):1095–1098. *(origin of S = 2.69 ± 0.19)*
- Slessarev et al. (2007). Prospective targeting and control of end-tidal CO₂ and O₂.
  *J Physiol* 581(3):1207–1219.
- West (2012). *Respiratory Physiology: The Essentials.*
- Lumb & Thomas (2020). *Nunn's Applied Respiratory Physiology.*
- Read (1967). A clinical method for assessing the ventilatory response to CO₂.
  *Australas Ann Med* 16(1):20–32.
- Bellville et al. (1979). Central and peripheral chemoreflex loop gain. *J Appl Physiol*
  46(4):843–853. *(central time constant)*
- Tansley et al. (1998). Human ventilatory response to 8 h of euoxic hypercapnia. *J Appl
  Physiol* 84(2):431–434. *(fast/central components)*
- Cunningham, Robbins & Wolff (1986). Integration of respiratory responses… *Handbook of
  Physiology, Sect. 3, Vol. II.* *(peripheral time constant)*

The **companion technical note** gives the full, cited derivation of both models (dead-space
basis, the two directions, the HCVR slope, time dependence and two-threshold structure, apnea
and anaesthesia, and validation against reported challenges).

## License

MIT — see [LICENSE](LICENSE).
