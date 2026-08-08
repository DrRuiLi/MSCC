# chemform_isotope_label

replace a specified number of atoms of an element with an isotope in a
chemical formula

## Usage

``` r
chemform_isotope_label(chemform = demo_chemform, ele = "[13]C", count = 1)
```

## Arguments

- chemform:

  chemical formula

- ele:

  isotope notation, e.g. `\\[13\\]C`, `\\[2\\]H`

- count:

  number of atoms to replace with the isotope

## Value

character, the labeled chemical formula

## Examples

``` r
chemform_isotope_label("C6H12O6", "[13]C", 3)
#> [1] "[13]C3C3H12O6"
```
