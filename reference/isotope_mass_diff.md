# isotope_mass_diff

calculate the mass difference of a given isotope to its major (most
abundant) form

## Usage

``` r
isotope_mass_diff(isotope)
```

## Arguments

- isotope:

  string, isotope notation such as `\\[13\\]C`, `\\[2\\]H`

## Value

numeric, mass difference (M\\13\\C - M\\12\\C)

## Examples

``` r
isotope_mass_diff("[13]C")
#> [1] 1.003355
isotope_mass_diff("[2]H")
#> [1] 1.006277
```
