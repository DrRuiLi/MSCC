# CFM_shiny

``` r

library(MSCC)
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
#> Loading required package: magrittr
#> 
#> 
#> Attaching package: 'magrittr'
#> 
#> 
#> The following object is masked from 'package:purrr':
#> 
#>     set_names
#> 
#> 
#> The following object is masked from 'package:tidyr':
#> 
#>     extract
```

## 1. Introduction

[`shiny_vis_cfm()`](https://drruili.github.io/MSCC/reference/shiny_vis_cfm.md)
is an interactive Shiny viewer for `CFM_data` objects (predicted /
annotated spectra from CFM-ID).

Layout:

| panel            | content                                         |
|------------------|-------------------------------------------------|
| Multi-CE spectra | plotly spectra (energy0/1/2 → CE 10/20/40)      |
| Precursor        | precursor molecule (`vis_smiles`)               |
| Fragment         | selected fragment molecule                      |
| Info             | precursor / fragment (id, SMILES, formula, m/z) |

Click a peak in the spectra to update Fragment and Info.

## 2. Usage

### 2.1 Create a CFM demo

`CFM_data` holds:

- `@peak_assignment` — peaks per CE with optional `fragment_id` /
  `fragment_score`
- `@fragment_define` — fragment SMILES and `fragment_mz`
- `@fragment_transition` — fragment graph edges (when from annotate)

Build one from a SMILES string via CFM-ID (Docker required):

``` r

cfm <- get_CFM_data_from_smiles(
  smiles = "NCC(O)=O",
  compound_id = "glycine",
  adduct = "[M+H]+",
  check_cache = TRUE,
  cache_dir = tempdir()
)
cfm
```

Small molecules such as glycine (`NCC(O)=O`) are useful for a quick
demo.

### 2.2 CFM shiny

Launch the viewer:

``` r

shiny_vis_cfm(cfm)
```

Optionally plot the multi-CE spectra alone:

``` r

plotly_CFM_spectra(cfm)
```
