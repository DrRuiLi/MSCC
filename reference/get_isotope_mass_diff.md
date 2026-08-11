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
#>      [13]C     [13]C2     [13]C3     [13]C4       [2]H      [18]O      [17]O 
#>  1.0033548  2.0067097  3.0100645  4.0134194  1.0062767  2.0042458  1.0042169 
#>     [18]O2 [17]O[18]O      [34]S      [33]S      [36]S 
#>  4.0084916  3.0084627  1.9957961  0.9993878  3.9950101 
get_isotope_mass_diff(element = c("C", "N", "S"), threshold = 0.01)
#>     [13]C     [15]N     [34]S     [33]S     [36]S 
#> 1.0033548 0.9970350 1.9957961 0.9993878 3.9950101 
```
