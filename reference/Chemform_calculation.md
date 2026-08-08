# Chemform Calculation

Chemform Calculation

Chemform Calculation

## Usage

``` r
chemform_sum(..., return = c("chemform", "matrix"))

chemform_multi(
  chemform = MSCC::demo_chemform,
  multi = 1,
  return = c("matrix", "chemform")
)
```

## Arguments

- ...:

  chemform list or vector

- chemform:

  chemform

- multi:

  number to multi

## Value

chemform

chemform

## Functions

- `chemform_sum()`: chemform_sum

- `chemform_multi()`: chemform_multi

## Examples

``` r
chemform_sum(demo_chemform)
#> [1] "C1242H1542Cl48N299O225S29Na5Br1P1F24Si1[2]H88[15]N5[13]C1"
chemform_multi("C6H12O6",1)
#>         C  H O
#> C6H12O6 6 12 6
```
