# Plotly multi-CE spectrum for CFM_data

Plotly multi-CE spectrum for CFM_data

## Usage

``` r
plotly_CFM_spectra(cfmd, source = "cfm_spectra")
```

## Arguments

- cfmd:

  A `CFM_data` object.

- source:

  Plotly event source id (default `"cfm_spectra"`).

## Value

A plotly object with clickable peaks (customdata = fragment_id).
