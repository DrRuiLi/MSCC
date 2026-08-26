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
  mzabs = NULL,
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

  Allowed deviation in ppm. `NULL` is the same as `Inf` (unused).
  Default `5`. Combined with `mzabs` by taking the **tighter** (smaller)
  window on the **neutral** mass used for MCP; they are not added. See
  **Mass tolerance**.

- mzabs:

  Allowed absolute deviation in Dalton. `NULL` is the same as `Inf`
  (unused). Default `NULL` (ppm-only). Combined with `ppm` by taking the
  tighter window (see **Mass tolerance**).

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

### Mass tolerance

Same rule as
[`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md):
`NULL` and `Inf` mean unused. `ppm` is converted to Dalton in R, then
the **minimum** of the two windows is passed to C++:

`abs_error = min(ppm * |M| * 1e-6, mzabs)`

on the converted **neutral** mass `M`. Defaults are `ppm = 5` and
`mzabs = NULL` (5 ppm). Set `ppm = NULL` (or `Inf`) for an absolute
window only.

## See also

[`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md),
[`chemform_mz()`](https://drruili.github.io/MSCC/reference/chemform_mz.md),
[`chemform_check_seven_golden_rules()`](https://drruili.github.io/MSCC/reference/chemform_check_seven_golden_rules.md)
