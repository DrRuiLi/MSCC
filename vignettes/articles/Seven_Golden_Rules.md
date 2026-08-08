# Seven Golden Rules (Kind & Fiehn) and `MassTools::calcMF`

Heuristic filters for molecular formulas from accurate mass, from:

> Kind T, Fiehn O (2007) *Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry.* BMC Bioinformatics 8:105.  
> https://doi.org/10.1186/1471-2105-8-105

## MSCC API

Use **`chemform_check_seven_golden_rules()`** to apply Rules **#1, #2 (Senior3), #4, #5, #6** to formula strings (no MassTools dependency; same thresholds as `MassTools::calcMF` post-filters). Rules **#3** (isotope) and **#7** (TMS) are not applied.

```r
chemform_check_seven_golden_rules(c("C6H12O6", "CH6N2"), return = "chemform")
chemform_check_seven_golden_rules(c("C6H12O6", "CH6N2"), return = "valid")
chemform_check_seven_golden_rules(c("C6H12O6", "CH6N2"), return = "all")

# Optional post-filter on MCP candidates:
chemform_decompose_mass(180.0634, check_rule = TRUE)
chemform_decompose_mz(181.0707, charge = 1, check_rule = TRUE)
```

| `return` | Result |
|----------|--------|
| `"chemform"` | Character vector of formulas that pass enabled rules (default) |
| `"valid"` | Logical vector, one per input |
| `"all"` | `data.frame` with input `chemform`, ratios/scores, per-rule flags, and `valid` |

`MassTools::calcMF` enumerates candidates with `Rdisop::decomposeMass`, then can apply **many** of these filters (defaults: all implemented golden-rule filters **on**). Rules **#3** (isotope pattern) and **#7** (TMS) are **not** applied inside `calcMF`.

---

## The seven rules (summary)

| # | Name | Purpose |
|---|------|---------|
| 1 | Element count limits | Cap atom counts by mass range (natural-product–like space) |
| 2 | LEWIS / SENIOR | Reject formulas that cannot form a chemically valid graph |
| 3 | Isotopic pattern | Score / reject by measured M+1, M+2 abundances |
| 4 | H/C ratio | Keep typical hydrogen/carbon ratios |
| 5 | Heteroatom / C ratios | Cap N/C, O/C, P/C, S/C (and related) |
| 6 | Multi-element probability | Reject improbable **combinations** of high N, O, P, S |
| 7 | TMS check | For GC–MS silylation: strip TMS, then re-apply rules |

Apply to **neutral** formulas (correct adducts / charge first).

---

### Rule 1 — Restrictions on element numbers

Hard upper bounds on C, H, N, O, P, S, … by mass window, derived from DNP / Wiley-type formula sets (not bare `mass / element_mass`).

**`calcMF` (`maxCounts`):** uses DNP-style caps (higher of DNP vs Wiley in the paper’s Table 1 for CHNOPS; halogens as in DNP):

| Mass \|m/z\| | C | H | N | O | P | S | F | Cl | Br |
|-------------|---|---|---|---|---|---|---|----|----|
| &lt; 500 | 29 | 72 | 10 | 18 | 4 | 7 | 15 | 8 | 5 |
| &lt; 1000 | 66 | 126 | 25 | 27 | 6 | 8 | 16 | 11 | 8 |
| &lt; 2000 | 115 | 236 | 32 | 63 | 6 | 8 | 16 | 11 | 8 |
| &lt; 3000 | 162 | 208 | 48 | 78 | 6 | 9 | 16 | 11 | 8 |
| ≥ 3000 | no maxCounts filter |

Also: `Filters$minElements` / `maxElements` (passed into `decomposeMass`).

---

### Rule 2 — LEWIS and SENIOR

**LEWIS (octet / valence):** atoms should have filled valence shells; strict enforcement rejects odd-electron species (e.g. many nitroso / radicals). Often relaxed via DBE / parity.

**SENIOR (molecular graph):** for a formula to exist as a connected graph,

1. Sum of valences (or count of odd-valence atoms) is even  
2. Sum of valences ≥ 2 × maximum valence  
3. Sum of valences ≥ 2 × (number of atoms − 1)

**`calcMF` (`SeniorRule`):** implements Senior’s **third** condition only:

```text
SENIOR3 = Σ (count × max_valence) − 2 × (n_atoms − 1)
```

Default keep if `SENIOR3 >= 0` (`Filters$SENIOR3 = 0`). Uses maximum valence states per element (as in Kind & Fiehn).

Extra related filters (not numbered golden rules): `Filters$DBErange` (default −5…40), `Filters$parity` (`"e"` / `"o"`).

---

### Rule 3 — Isotopic pattern filter

Compare predicted vs measured isotope abundances (especially M+1, M+2). Orthogonal to mass error; essential above ~200 Da for CHNOPS uniqueness.

**Not in `calcMF`.** Use isotope-aware tools separately (e.g. `Rdisop::decomposeIsotopes`, enviPat / MSCC isotope helpers).

---

### Rule 4 — Hydrogen / carbon ratio

Common range (≈99.7% of Wiley formulas): **0.2 &lt; H/C &lt; 3.1**.  
Extended (≈99.99%): 0.1–6. Extremes (fullerenes, etc.) need the rule disabled.

**`calcMF` (`HCratio`):** keeps `0.2 < H/C < 3.1`. Skipped if H or C absent (`NA`).

---

### Rule 5 — Heteroatom / carbon ratios (NOPS)

Common upper bounds (Wiley, Table 2 of Kind & Fiehn):

| Ratio | Common max |
|-------|------------|
| N/C | &lt; 1.3 |
| O/C | &lt; 1.2 |
| P/C | &lt; 0.3 |
| S/C | &lt; 0.8 |

(Paper also lists F/C, Cl/C, Br/C, Si/C; `calcMF` only checks N, O, P, S vs C.)

**`calcMF` (`moreRatios`):** same four caps. Ignored if no carbon.

---

### Rule 6 — Element probability (multiple high HNOPS)

Rejects formulas with **several** heteroatoms simultaneously high (e.g. `C26H28N17O1P3S8` can pass #4–5 but fail #6).

Paper Table 3 (compounds &lt; 2000 Da), mirrored in **`calcMF` (`elementHeuristic`)**:

| Condition | Reject if |
|-----------|-----------|
| N,O,P,S all &gt; 1 | N ≥ 10 or O ≥ 20 or P ≥ 4 or S ≥ 3 |
| N,O,P all &gt; 3 | N ≥ 11 or O ≥ 22 or P ≥ 7* |
| O,P,S all &gt; 1 | O ≥ 14 or P ≥ 3 or S ≥ 3 |
| N,P,S all &gt; 1 | N ≥ 4 or P ≥ 3 or S ≥ 3 |
| N,O,S all &gt; 6 | N ≥ 19 or O ≥ 14 or S ≥ 8 |

\*Paper Table 3 lists P &lt; 6 for NOP; MassTools uses `P >= 7` as the fail threshold.

If fewer than three of N/O/P/S exceed count 1, the heuristic passes.

---

### Rule 7 — TMS check (GC–MS)

After silylation, subtract TMS units (`C3H8Si` per group), then re-apply ratio / valence rules on the **native** formula. Si isotope pattern helps count TMS groups.

**Not in `calcMF`.** Handle derivatization (like adducts) before calling the calculator.

---

## What `MassTools::calcMF` turns on by default

```r
calcMF(
  mz, z = 1, ppm = 5,
  maxCounts = TRUE,       # Rule 1
  SeniorRule = TRUE,      # Rule 2 (Senior #3)
  HCratio = TRUE,         # Rule 4
  moreRatios = TRUE,      # Rule 5
  elementHeuristic = TRUE,# Rule 6
  Filters = list(
    DBErange = c(-5, 40),
    minElements = "C0",
    maxElements = "C99999",
    parity = "e",
    maxCounts = TRUE,
    SENIOR3 = 0,
    HCratio = TRUE,
    moreRatios = TRUE,
    elementHeuristic = TRUE
  )
)
```

| Golden rule | In `calcMF`? | Argument / note |
|-------------|--------------|-----------------|
| #1 Element counts | Yes | `maxCounts` |
| #2 LEWIS / SENIOR | Partial | `SeniorRule` / `SENIOR3`; plus DBE & parity |
| #3 Isotope pattern | No | — |
| #4 H/C | Yes | `HCratio` → 0.2–3.1 |
| #5 NOPS/C | Yes | `moreRatios` |
| #6 HNOPS probability | Yes | `elementHeuristic` |
| #7 TMS | No | Strip TMS upstream |

Set a filter flag to `FALSE` (or omit the corresponding `Filters` entry / set it `NULL` where documented) to disable that step.

---

## Practical notes

1. Neutralize adducts (`[M+H]+` → subtract H, etc.) before trusting valence / ratio filters.  
2. Defaults target **small-molecule / natural-product–like** space; peptides and exotic chemistries may need looser or disabled rules.  
3. For LC–MS formula ID, combine #1–2 and #4–6 with **isotope fit (#3)** and database lookup when available.  
4. Open-access article and lab notes: [BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-105), [Fiehn Lab – Seven Golden Rules](https://www.fiehnlab.ucdavis.edu/projects/seven-golden-rules).
