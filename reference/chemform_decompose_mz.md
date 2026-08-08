# MCP mass decomposition from ion m/z

Convert ion m/z + charge to neutral exact mass (same electron-mass
convention as
[`chemform_mz()`](https://drruili.github.io/MSCC/reference/chemform_mz.md)),
then call
[`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md).

## Usage

``` r
chemform_decompose_mz(
  mz,
  charge = 0,
  ppm = 5,
  mzabs = 1e-04,
  elements = c("C", "H", "N", "O", "P", "S"),
  min_elements = NULL,
  max_elements = NULL,
  check_rule = FALSE
)
```

## Arguments

- mz:

  Ion m/z (scalar or vector). If `charge = 0`, treated as neutral mass.

- charge:

  Integer charge `z` (e.g. `+1`, `+2`, `-1`). Use `0` for neutral input.

- ppm:

  Allowed deviation in ppm (applied on the neutral mass used for MCP).

- mzabs:

  Allowed absolute deviation in Dalton.

- elements:

  Character vector of allowed elements.

- min_elements:

  See
  [`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md).

- max_elements:

  See
  [`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md).

- check_rule:

  If `TRUE`, keep only candidates that pass
  [`chemform_check_seven_golden_rules()`](https://drruili.github.io/MSCC/reference/chemform_check_seven_golden_rules.md)
  (passed through to
  [`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md)).
  Default `FALSE`. Disabled automatically when any `min_elements` count
  is negative (see
  [`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md)).

## Value

A data.frame with columns `formula`, `exactmass`, `mz`, `ppm`, `charge`,
and `mz_target`. Rows are sorted by increasing `abs(ppm)` vs the input
m/z.

## Details

Neutral mass conversion when `charge != 0`:
`M = mz * abs(charge) + e * charge`, with `e = 0.00054857990943`.

## See also

[`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md),
[`chemform_mz()`](https://drruili.github.io/MSCC/reference/chemform_mz.md),
[`chemform_check_seven_golden_rules()`](https://drruili.github.io/MSCC/reference/chemform_check_seven_golden_rules.md)
