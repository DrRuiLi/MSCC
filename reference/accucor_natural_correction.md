# AccuCor-style natural isotope correction for ratio matrix

AccuCor-style natural isotope correction for ratio matrix

## Usage

``` r
accucor_natural_correction(
  raw_ratio,
  formula,
  iso_ele = "[13]C",
  purity = 1,
  Resolution = 140000,
  ResDefAt = 200
)
```

## Arguments

- raw_ratio:

  Numeric matrix (rows = isotopologues, cols = samples).

- formula:

  Molecular formula string.

- iso_ele:

  Tracer isotope, e.g. `"[13]C"`, `"[2]H"`.

- purity:

  Tracer purity.

- Resolution:

  MS resolving power.

- ResDefAt:

  m/z where resolution is defined.

## Value

Numeric matrix with corrected ratios.
