# Convert CFM Data to Spectra Object

Converts CFM (Competitive Fragmentation Modeling) peak assignment data
to a `Spectra` object with spectra at three collision energy levels (10,
20, 40 eV for energy0/1/2).

## Usage

``` r
get_CFM_data_Spectra(cfmd)
```

## Arguments

- cfmd:

  A `CFM_data` object containing `peak_assignment`. Typically from
  [`read_CFM_predict_result`](https://drruili.github.io/MSCC/reference/CFM.md)
  /
  [`CFM_annotate_by_predict`](https://drruili.github.io/MSCC/reference/CFM.md).

## Value

A `Spectra` object containing MS/MS spectra at three CE levels.

## Details

Convert CFM_data peak assignments to a Spectra object
