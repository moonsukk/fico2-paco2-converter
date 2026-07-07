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

You pick a direction and are prompted for the inputs:

- **Direction 1 (PaCO₂ → FiCO₂):** enter target PaCO₂ and baseline PaCO₂.
- **Direction 2 (FiCO₂ → PaCO₂):** enter FiCO₂; you're asked whether you know the baseline
  PaCO₂. If not, the resting fallback (`VA_base = 4.2`, `PaCO2_base = 40`) is used.

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

### Jupyter notebook

`FiCO2_PaCO2_conversion.ipynb` walks through the derivation, both directions, a reference
table, and a sensitivity analysis over the HCVR slope `S`.

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
- HCVR treated as linear over the working range (≈ 40–80 mmHg).
- A single population-mean slope `S = 2.69`; it will not match any individual — see the
  sensitivity analysis in the notebook.
- Published `S` are usually minute-ventilation slopes; used here for the alveolar slope,
  valid while dead-space ventilation changes little.

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

See `review_derivation` (simplified, fixed-slope form) and `technical_note` (reversible,
multi-unit, adjustable slope) for the full derivation.

## License

MIT — see [LICENSE](LICENSE).
