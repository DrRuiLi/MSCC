# chemform_adduct

get chemical formula with adduct

## Usage

``` r
chemform_adduct(
  chemform = MSCC::demo_chemform,
  adduct = "[M+H]+",
  value = c("mz", "chemform", "all")
)
```

## Arguments

- chemform:

  Chemical formula.

- adduct:

  Adduct form, such as "\[M+H\]+".

- value:

  Output type: m/z, formula, or full result table.

## Value

mz or df
