# chemform_calc

chemform calculator, calc "+" or "-" for vector operation, which return
a vetor with length of m. calc ".+" or ".-" for matrix operation, which
return a matrix m x n

## Usage

``` r
chemform_calc(
  chemform1 = MSCC::demo_chemform,
  chemform2 = rev(MSCC::demo_chemform),
  calc = "+",
  return = c("matrix", "chemform")
)
```

## Arguments

- chemform1:

  vector of chemform, length as m

- chemform2:

  vector of chemform, length as n

- calc:

  - or -

## Value

chemform vector or matrix
