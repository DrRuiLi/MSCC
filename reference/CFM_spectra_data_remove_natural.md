# Remove Natural Isotope Contribution from Combined CFM Spectra

Subtracts estimated natural isotopologue intensity from combined
spectrum peaks (`sp.id = "combined_sp"`) using an isoform map and a
natural-abundance scaling factor. Peaks with non-positive intensity
after subtraction are set to zero; fragment groups that become all-zero
are dropped from the combined spectrum.

## Usage

``` r
CFM_spectra_data_remove_natural(sp.data, natural.ratio, if.map)
```

## Arguments

- sp.data:

  Data frame of spectrum peaks that may include a combined spectrum from
  [`CFM_spectra_data_int_weight`](https://drruili.github.io/MSCC/reference/CFM_spectra_data_int_weight.md)
  or
  [`CFM_spectra_data_merge`](https://drruili.github.io/MSCC/reference/CFM_spectra_data_merge.md).

- natural.ratio:

  Numeric scaling factor for the natural isotope contribution.

- if.map:

  Object with an `isoform.map` slot used to derive the natural
  isotopologue distribution.

## Value

A data frame with natural contribution removed from combined peaks.

## Details

Remove Natural Isotope Contribution from Combined CFM Spectra
