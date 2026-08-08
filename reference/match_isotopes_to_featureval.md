# Append feature intensity matrix to isotope matches

Append feature intensity matrix to isotope matches

## Usage

``` r
match_isotopes_to_featureval(matched_table, featureval)
```

## Arguments

- matched_table:

  Output from
  [`match_isotopes_to_featuredef()`](https://drruili.github.io/MSCC/reference/match_isotopes_to_featuredef.md).

- featureval:

  Numeric matrix from
  [`xcms::featureValues()`](https://rdrr.io/pkg/xcms/man/XCMSnExp-peak-grouping-results.html).

## Value

Data frame with matched rows and joined intensity columns.
