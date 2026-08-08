# Weight CFM Spectra Intensities by Fragment Group

Computes intensity-weighted mean isotopologue intensities for each
fragment group across spectra, then appends a combined spectrum
(`sp.id = "combined_sp"`) to the input data.

## Usage

``` r
CFM_spectra_data_int_weight(sp.data, iso_count)
```

## Arguments

- sp.data:

  Data frame of annotated spectrum peaks, typically from
  [`CFM_annotate_isotopologues`](https://drruili.github.io/MSCC/reference/MSIP.md),
  with columns such as `fragment_group`, `sp.id`, `iso`, `intensity`,
  and `mz`.

- iso_count:

  Maximum isotope count (Mn) to include, as an integer.

## Value

A data frame with the original peaks plus weighted combined peaks.

## Details

Weight CFM Spectra Intensities by Fragment Group
