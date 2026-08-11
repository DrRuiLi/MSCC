# MSCC

**MSCC** (Mass Spectrometry Chemical Calculation) provides
chemical-formula arithmetic, isotope and adduct handling, molecular
graph helpers, and CFM-ID spectrum prediction/annotation for
metabolomics and mass spectrometry workflows.

📚 **Documentation:** <https://drruili.github.io/MSCC/>

## Installation

Install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("WallFcerLR/MSCC")
```

Or with `remotes`:

``` r

# install.packages("remotes")
remotes::install_github("WallFcerLR/MSCC")
```

## Quick start

``` r

library(MSCC)

# Parse formulas and compute exact mass / m/z
chemform_parse("C6H12O6")
chemform_mz("C6H12O6")
chemform_mz("C6H12O6", charge = 1)

# Adducts and isotope labels
chemform_adduct("C6H12O6", "[M+H]+")
chemform_isotope_label("C6H12O6", "[13]C", 3)

# Formula candidates from accurate mass or m/z
chemform_decompose_mass(180.0634, ppm = 5)
chemform_decompose_mz(181.0707, charge = 1, ppm = 5)

# Heuristic filters (Seven Golden Rules)
chemform_check_seven_golden_rules(c("C6H12O6", "CH6N2"))
```

## Features

| Area | Highlights |
|----|----|
| **Formulas** | Parse, format, sum/multiply, exact mass and charged m/z |
| **Adducts / isotopes** | Adduct tables, isotope mass diffs, labeled formulas, enviPat patterns |
| **Decomposition** | MCP-based [`chemform_decompose_mass()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mass.md) / [`chemform_decompose_mz()`](https://drruili.github.io/MSCC/reference/chemform_decompose_mz.md) with optional golden-rule filtering |
| **Molecules** | `Molecule_igraph` (SDF + igraph + isotopomers), SMILES/SDF bridges, visualization |
| **CFM-ID** | Docker-backed prediction/annotation (`CFM_data`), Spectra I/O, Shiny viewer |
| **Optional chem** | RDKit / RXNMapper via reticulate; natural-abundance correction (`accucor`) |

## Conventions

- **Isotope notation:** `[13]C`, `[2]H`, `[15]N`, `[18]O`, `[34]S`
- **Adduct strings:** prefer `[M+H]+` / `[M-H]-`
- **Polarity:** `0` = negative, `1` = positive
- **CFM collision energies:** `energy0` → CE 10, `energy1` → CE 20,
  `energy2` → CE 40

## Documentation

Full reference and articles are built with
[pkgdown](https://pkgdown.r-lib.org/) and published to GitHub Pages:

- Site: <https://drruili.github.io/MSCC/>
- Source: <https://github.com/WallFcerLR/MSCC>
- Issues: <https://github.com/WallFcerLR/MSCC/issues>

Rebuild the site locally:

``` r

pkgdown::build_site()
```

## License

GPL-2 © [Rui Li](mailto:rli@sinh.ac.cn)
