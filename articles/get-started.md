# Get started

``` r

library(MSCC)
```

## Installation

Install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("drruili/MSCC")
```

Or with `remotes`:

``` r

# install.packages("remotes")
remotes::install_github("drruili/MSCC")
```

## Formulas, mass, and m/z

Parse a chemical formula to an element matrix, then compute exact mass
or charged m/z:

``` r

chemform_parse("C6H12O6")
#>         C  H O
#> C6H12O6 6 12 6
chemform_mz("C6H12O6")
#> [1] 180.0634
chemform_mz("C6H12O6", charge = 1)
#> [1] 180.0628
```

Sum or multiply formulas when building larger species:

``` r

chemform_sum("C6H12O6", "H2O")
#> [1] "C6H14O7"
chemform_multi("CH2", 3)
#>      C H
#> C1H2 3 6
```

## Adducts and isotopes

Prefer adduct strings like `[M+H]+` / `[M-H]-`. Isotope labels use
bracket notation (`[13]C`, `[2]H`, …):

``` r

chemform_adduct("C6H12O6", "[M+H]+", value = "chemform")
#> [1] "C6H13O6"
chemform_isotope_label("C6H12O6", "[13]C", 3)
#> [1] "[13]C3C3H12O6"
```

Optional enviPat isotope patterns (requires the **enviPat** package):

``` r

chemform_isotopes_pattern_enviPat("C6H12O6")
#> # A tibble: 6 × 4
#> # Rowwise: 
#>   formula         m.z abundance isotope_element
#>   <chr>         <dbl>     <dbl> <chr>          
#> 1 C6H12O6        180.   100     ""             
#> 2 C5H12O6[13]C1  181.     6.49  "[13]C1"       
#> 3 C6H12O5[18]O1  182.     1.23  "[18]O1"       
#> 4 C6H12O5[17]O1  181.     0.229 "[17]O1"       
#> 5 C4H12O6[13]C2  182.     0.175 "[13]C2"       
#> 6 C6H11O6[2]H1   181.     0.138 "[2]H1"
```

## Formula candidates from accurate mass

Decompose a neutral exact mass or an ion m/z into candidate formulas.
Optionally apply Seven Golden Rules filters:

``` r

chemform_decompose_mass(180.0634, ppm = 5)
#>     formula exactmass         ppm mass_target
#> 1   C6H12O6  180.0634 -0.06606562    180.0634
#> 2   C5H6N7O  180.0634 -0.09506652    180.0634
#> 3 CH17N4PS2  180.0632 -0.97446788    180.0634
#> 4  H23O2P3S  180.0632 -1.33188644    180.0634
#> 5   CH9N8OP  180.0637  1.62985926    180.0634
#> 6 C2H15NO6P  180.0637  1.65886016    180.0634
#> 7 C5H14N3S2  180.0629 -2.69939366    180.0634
#> 8  C7H16OS2  180.0643  4.75716886    180.0634
chemform_decompose_mz(181.0707, charge = 1, ppm = 5)
#>     formula exactmass       mz        ppm charge mz_target
#> 1   C6H13O6  181.0712 181.0707 -0.1957462      1  181.0707
#> 2   C5H7N7O  181.0712 181.0707 -0.2245858      1  181.0707
#> 3 CH18N4PS2  181.0710 181.0705 -1.0990950      1  181.0707
#> 4  H24O2P3S  181.0710 181.0704 -1.4545253      1  181.0707
#> 5  CH10N8OP  181.0715 181.0710  1.4907442      1  181.0707
#> 6 C2H16NO6P  181.0715 181.0710  1.5195837      1  181.0707
#> 7 C5H15N3S2  181.0707 181.0702 -2.8144250      1  181.0707
#> 8  C7H17OS2  181.0721 181.0715  4.6006565      1  181.0707

chemform_check_seven_golden_rules(c("C6H12O6", "CH6N2"))
#> [1] "C6H12O6"
```

See the [Seven Golden
Rules](https://drruili.github.io/MSCC/articles/articles/Seven_Golden_Rules.md)
article for rule details.

## Molecules from SMILES

Build a `Molecule_igraph` (SDF + igraph + optional isotopomers). Default
demo molecule is glycine:

``` r

mol <- get_Molecule_igraph_from_smiles("NCC(O)=O", id = "glycine")
mol@molecule_info
```

## CFM-ID spectra (optional)

Prediction and annotation need a working CFM-ID Docker setup. Small
molecules such as glycine are useful for a quick demo:

``` r

cfm <- get_CFM_data_from_smiles(
  smiles = "NCC(O)=O",
  compound_id = "glycine",
  adduct = "[M+H]+",
  check_cache = TRUE,
  cache_dir = tempdir()
)
shiny_vis_cfm(cfm)
```

More detail: the
[CFM_shiny](https://drruili.github.io/MSCC/articles/CFM_shiny.md)
vignette.

## Conventions

| Topic            | Convention                                              |
|------------------|---------------------------------------------------------|
| Isotope notation | `[13]C`, `[2]H`, `[15]N`, `[18]O`, `[34]S`              |
| Adduct strings   | Prefer `[M+H]+` / `[M-H]-`                              |
| Polarity         | `0` = negative, `1` = positive                          |
| CFM energies     | `energy0` → CE 10, `energy1` → CE 20, `energy2` → CE 40 |

## Next steps

- Reference index: function-level help on the package site
- Articles: formula decomposition, Seven Golden Rules, CFM viewer, RDKit
  helpers
