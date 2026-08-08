# Merge CFM Spectra Isotopologue Ratios by Fragment Group

Merges per-spectrum isotopologue ratios within each fragment group using
intensity-weighted means, and appends combined peaks
(`sp.id = "combined_sp"`, `merged = TRUE`) with summary metrics
(`int_sum`, `peaks_count`, `icc`, `cos`).

## Usage

``` r
CFM_spectra_data_merge(sp.data, iso_count)
```

## Arguments

- sp.data:

  Data frame of annotated spectrum peaks with fragment-group ratios,
  typically from
  [`CFM_annotate_isotopologues`](https://drruili.github.io/MSCC/reference/MSIP.md),
  including columns such as `fragment_group`, `sp.id`, `iso_count`,
  `ratio`, `int_sum`, and `mz`.

- iso_count:

  Maximum isotope count (Mn) to include, as an integer.

## Value

A data frame with original peaks plus merged combined peaks and metrics.

## Details

Merge CFM Spectra Isotopologue Ratios by Fragment Group
