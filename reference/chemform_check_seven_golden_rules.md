# Check formulas against Kind & Fiehn seven golden rules (#1, \#2, \#4–#6)

Reimplements the MassTools::calcMF post-filters (no MassTools
dependency). Rules \#3 (isotope pattern) and \#7 (TMS) are not applied.

## Usage

``` r
chemform_check_seven_golden_rules(
  chemform,
  mass = NULL,
  maxCounts = TRUE,
  SeniorRule = TRUE,
  HCratio = TRUE,
  moreRatios = TRUE,
  elementHeuristic = TRUE,
  DBErange = NULL,
  parity = NULL,
  SENIOR3_min = 0,
  return = c("chemform", "valid", "all")
)
```

## Arguments

- chemform:

  Character vector of chemical formulas (neutral preferred).

- mass:

  Exact mass(es) for Rule \#1 mass bins. Scalar recycled, or length
  equal to `chemform`. If `NULL`, monoisotopic mass is computed from the
  formula.

- maxCounts:

  Logical; apply element count caps by mass range (Rule \#1).

- SeniorRule:

  Logical; apply Senior's third theorem (Rule \#2).

- HCratio:

  Logical; apply H/C common range 0.2–3.1 (Rule \#4).

- moreRatios:

  Logical; apply N/O/P/S vs C ratio caps (Rule \#5).

- elementHeuristic:

  Logical; apply multi-HNOPS probability check (Rule \#6).

- DBErange:

  Optional numeric length-2 `c(min, max)` for DBE filter; `NULL`
  disables.

- parity:

  Optional `"e"` or `"o"` (integer vs half-integer DBE); `NULL`
  disables.

- SENIOR3_min:

  Minimum Senior3 score to pass (default `0`).

- return:

  One of `"chemform"` (passing formulas), `"valid"` (logical vector), or
  `"all"` (diagnostic data.frame).

## Value

Depends on `return`:

- `"chemform"`: character vector of formulas that pass enabled rules

- `"valid"`: logical vector, one per input

- `"all"`: data.frame with input `chemform`, ratios/scores, per-rule
  flags, and `valid`

## References

Kind T, Fiehn O (2007) Seven Golden Rules for heuristic filtering of
molecular formulas obtained by accurate mass spectrometry. BMC
Bioinformatics 8:105.

## See also

[`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md),
[`chemform_parse()`](https://drruili.github.io/MSCC/reference/chemform_parse.md)
