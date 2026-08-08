# MCP-only mass decomposition (neutral exact mass)

Enumerate candidate sum formulas that match a **neutral exact mass**
within an error window. Uses an imslib/Rdisop-style Money Changing
Problem (MCP) solver and intentionally skips isotope distribution
calculation/ranking.

## Usage

``` r
chemform_decompose_mass(
  mass,
  ppm = 5,
  mzabs = 1e-04,
  elements = c("C", "H", "N", "O", "P", "S"),
  min_elements = NULL,
  max_elements = NULL,
  check_rule = FALSE
)
```

## Arguments

- mass:

  Neutral exact mass (scalar or vector).

- ppm:

  Allowed deviation in ppm.

- mzabs:

  Allowed absolute deviation in Dalton.

- elements:

  Character vector of allowed elements, e.g.
  `c("C","H","N","O","P","S")`.

- min_elements:

  Minimum element counts. Accepts `NULL` (defaults to 0 for each
  element), a named integer vector (names are element symbols), or a
  single formula string like `"C0H0N0"`. Counts may be **negative**
  (e.g. `c(H = -3)` or `"H-3"`) for replacement / difference formulas
  such as COOH → COONa (`H-1Na`). Internally the target mass is shifted
  so the MCP solver still enumerates non-negative relative counts.

- max_elements:

  Maximum element counts. Accepts `NULL` (defaults to 999999 for each
  element), a named integer vector (names are element symbols), or a
  single formula string like `"C999H999"`.

- check_rule:

  If `TRUE`, keep only candidates that pass
  [`chemform_check_seven_golden_rules()`](https://drruili.github.io/MSCC/reference/chemform_check_seven_golden_rules.md)
  (Rules \#1, \#2, \#4–#6). Default `FALSE`. Automatically disabled
  (with a message) when any `min_elements` count is negative, because
  the golden rules assume molecular formulas, not signed replacements.

## Value

A data.frame with columns `formula`, `exactmass`, `ppm`, and
`mass_target`. Rows are sorted by increasing `abs(ppm)`.

## Details

For ion m/z input with charge, use
[`chemform_decompose_mz()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mz.md).

## See also

[`chemform_decompose_mz()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mz.md),
[`chemform_check_seven_golden_rules()`](https://drruili.github.io/MSCC/reference/chemform_check_seven_golden_rules.md)

## Examples

``` r
# Neutral molecule (non-negative counts)
chemform_decompose_mass(180.0634, ppm = 5, check_rule = FALSE)
#>     formula exactmass         ppm mass_target
#> 1   C6H12O6  180.0634 -0.06606562    180.0634
#> 2   C5H6N7O  180.0634 -0.09506652    180.0634
#> 3 CH17N4PS2  180.0632 -0.97446788    180.0634
#> 4  H23O2P3S  180.0632 -1.33188644    180.0634
#> 5   CH9N8OP  180.0637  1.62985926    180.0634
#> 6 C2H15NO6P  180.0637  1.65886016    180.0634
#> 7 C5H14N3S2  180.0629 -2.69939366    180.0634
#> 8  C7H16OS2  180.0643  4.75716886    180.0634

# Replacement delta: COOH -> COONa is H-1Na (mass change Na - H)
dM <- chemform_mz("Na") - chemform_mz("H")
chemform_decompose_mass(
  mass = dM,
  elements = c("H", "C", "O", "Na"),
  min_elements = c(H = -3),
  max_elements = c(H = 10, Na = 3),
  check_rule = FALSE
)
#>   formula exactmass ppm mass_target
#> 1   H-1Na  21.98194   0    21.98194
```
