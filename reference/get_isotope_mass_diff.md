# get_isotope_mass_diff

For each formula / element fragment in `element`, runs
[`chemform_isotopes_pattern_enviPat()`](https://drruili.github.io/MSCC/reference/chemform_isotopes_pattern_enviPat.md)
and returns a named numeric vector of isotopologue mass differences vs
the monoisotopic peak. Names use compact isotope notation (e.g. \[13\]C,
\[13\]C2, \[34\]S). Accepts plain element symbols (`"C"`, `"S"`) or
counted fragments (`"C10"`).

## Usage

``` r
get_isotope_mass_diff(
  element = c("C10", "H10", "O5", "N5", "P3", "S3", "K", "Cl", "Br"),
  threshold = 1e-04
)
```

## Arguments

- element:

  Character vector of formulas / element symbols, e.g.
  `c("C10", "H10", "O2", "S")`.

- threshold:

  enviPat abundance cutoff passed as `thresh` (percent of the
  monoisotopic peak; default `0.0001`).

## Value

Named numeric vector: names are isotope labels, values are `mass_diff`.

## Examples

``` r
get_isotope_mass_diff(c("C10", "H10", "O2", "S"))
#> Warning: Element not found in elem_table: C10
#> Warning: Element not found in elem_table: H10
#> Warning: Element not found in elem_table: O2
#>     [33]S     [34]S 
#> 0.9993878 1.9957961 
get_isotope_mass_diff(element = c("C", "N", "S"), threshold = 0.01)
#> Error in get_isotope_mass_diff(element = c("C", "N", "S"), threshold = 0.01): All inputs must be named with allowed elements: C, H, O, N, P, S
```
