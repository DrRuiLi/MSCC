# Shiny app to visualize CFM_data spectra and fragments

Layout (left \| right):


    | p1 spectra | p2 precursor molecule |
    | p1 spectra | p3 fragment molecule  |
    | p1 spectra | p4 smiles / formula / mz |

## Usage

``` r
shiny_vis_cfm(cfmd, launch.browser = TRUE, host = "127.0.0.1")
```

## Arguments

- cfmd:

  A `CFM_data` object, typically from
  [`get_CFM_data_from_smiles`](https://drruili.github.io/MSCC/reference/get_CFM_data_from_smiles.md)
  or
  [`CFM_annotate_by_predict`](https://drruili.github.io/MSCC/reference/CFM.md).

- launch.browser:

  Logical; forward to
  [`shiny::shinyApp`](https://rdrr.io/pkg/shiny/man/shinyApp.html)
  options (default `TRUE`).

- host:

  Host binding for the app (default `"127.0.0.1"`).

## Value

A `shiny.appobj` (invisibly when run interactively via the app).

## Details

Click a peak in the multi-CE spectrum (p1) to update the fragment
molecule (p3) and the info panel (p4).

## Examples

``` r
if (FALSE) { # \dontrun{
cfm <- get_CFM_data_from_smiles("NCC(O)=O", compound_id = "glycine")
shiny_vis_cfm(cfm)
} # }
```
