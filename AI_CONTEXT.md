# AI_CONTEXT — MSCC

> Architectural map for AI assistants. Read this before exploring the
> codebase. Do **not** auto-update this file; propose changes and wait
> for explicit approval.

## Package identity

| Field   | Value                                                   |
|---------|---------------------------------------------------------|
| Name    | **MSCC** (Mass Spectrometry Chemical Calculation)       |
| Version | 1.1.0                                                   |
| Author  | Rui Li (`rli@sinh.ac.cn`)                               |
| Repo    | <https://github.com/WallFcerLR/MSCC>                    |
| Type    | R package (S4 + tidyverse + ChemmineR/igraph + Spectra) |

**Purpose:** Chemical-formula arithmetic, isotope/adduct handling for
MS, molecular graph representation, CFM-ID spectrum
prediction/annotation, and related metabolomics helpers (atom mapping,
natural-abundance correction, interactive viewers).

## Mental model

    SMILES / formula / SDF
            │
            ├─► chemform_*  ──► m/z, adducts, isotope patterns
            │
            ├─► Molecule_igraph (SDF + igraph + isotopomers)
            │         │
            │         ├─► atom maps / RXNMapper transfer
            │         └─► vis_smiles / visNetwork
            │
            └─► CFM_data (Docker CFM-ID)
                      │
                      ├─► Spectra / peak lists
                      └─► shiny_vis_cfm / plotly_CFM_spectra

## Core S4 classes

### `Molecule_igraph` (`R/Molecule_igraph-class.R`)

Slots: `molecule_info`, `sdf` (ChemmineR `SDF`), `igraph`, `isotopomer`
(data.frame).

- Build via `get_Molecule_igraph_from_smiles()` /
  `get_Molecule_igraph_from_sdf()`.
- Graph attrs via generics `vdata` / `edata` (defined in
  `R/00_graph_primitives.R`).
- Isotopomer helpers:
  [`Molecule_igraph_add_isotopomer()`](https://drruili.github.io/MSCC/reference/Molecule_igraph_add_isotopomer.md),
  atom transfer via MCS or RXNMapper.

### `CFM_data` (`R/dev_CFM.R`)

Slots:

| Slot | Role |
|----|----|
| `peak_assignment` | Peaks per CE (`energy0/1/2` → CE 10/20/40), optional `fragment_id` / `fragment_score` |
| `fragment_define` | Fragment SMILES, `fragment_mz`, polarity |
| `fragment_transition` | Fragment graph edges (annotate path) |

Typical entry:
[`get_CFM_data_from_smiles()`](https://drruili.github.io/MSCC/reference/get_CFM_data_from_smiles.md)
→ `CFM_predict` / `CFM_annotate*` (Docker). Viewer:
[`shiny_vis_cfm()`](https://drruili.github.io/MSCC/reference/shiny_vis_cfm.md).

## Module map (`R/`)

| File | Responsibility |
|----|----|
| `Chemform-function.R` | Parse/format formulas; `chemform_parse`, `chemform_get_ele`, `chemform_formate`, `chemform_mz`, `chemform_calc` |
| `Chemform_calc.R` | Formula sum/multiply (`chemform_sum`, `chemform_multi`) |
| `Adduct-function.R` | Adduct table logic; `chemform_adduct`, `chemform_adduct_check` |
| `isotope.R` | Isotope mass diffs, formula labeling (`chemform_isotope_label`) |
| `dev_envipat_lc8.R` | enviPat isotope patterns + lc8 m/z bridge |
| `migration_helpers.R` | Self-contained helpers (no MSdev); elem table, xcms isotope matching, polarity |
| `accucor_wrapper.R` | Natural-isotope correction (`accucor_natural_correction`) |
| `00_graph_primitives.R` | `vdata`/`edata`/`atom`/`get_element` generics; MCS/RXNMapper hooks |
| `Molecule_igraph-class.R` | S4 class + methods |
| `Molecule_igraph-functions.R` | Build/filter/vis isotopomers; atom-transfer matrices |
| `Molecule_vis.R` | SDF↔︎igraph, `vis_sdf_igraph` |
| `dev_Chemmine.R` | SMILES↔︎SDF, MCS atom maps, `vis_smiles`, isotopologue formatting |
| `dev_CFM.R` | CFM-ID Docker API, `CFM_data` I/O, Spectra bridges |
| `shiny_vis_cfm.R` | Shiny + plotly CFM viewer |
| `dev_rdkit.R` | reticulate RDKit (`get_RDKit_Chem`, formula/mass/SMARTS) |
| `dev_RXNmapper.R` | Reaction atom mapping via RXNMapper |
| `Bechmark.R` / `dev_NetID.R` | Dev/benchmark scratch (not primary API) |

## Packaged data (`data/`)

| Object                  | Use                                          |
|-------------------------|----------------------------------------------|
| `elem_table`            | Element/isotope masses & abundances          |
| `adduct.table`          | Adduct definitions + synonyms (`Adduct_Syn`) |
| `chem_formula_template` | Example / benchmark formulas                 |
| `isotopes_from_envipat` | enviPat isotope table                        |
| `smiles_map`            | SMILES lookup helper                         |

Access element table via
[`get_elem_table()`](https://drruili.github.io/MSCC/reference/get_elem_table.md)
when schema harmonization is needed.

## Conventions agents must respect

1.  **Isotope notation:** `[13]C`, `[2]H`, `[15]N`, `[18]O`, `[34]S`
    (not `13C` alone in chemform APIs).
2.  **Adduct strings:** Prefer `[M+H]+` / `[M-H]-`; many APIs also
    accept `+`/`-` or polarity `0`/`1`.
3.  **Polarity:** `0` = negative, `1` = positive
    (`get_polarity_from_adduct`, `get_polarity_suffix`).
4.  **CFM collision energies:** `energy0`→CE 10, `energy1`→CE 20,
    `energy2`→CE 40.
5.  **RDKit:** Import submodules explicitly
    ([`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md));
    do not invent a catch-all `get_RDKit()`.
6.  **CFM runtime:** Prediction/annotation needs a working CFM-ID Docker
    setup; local cache dirs are supported.
7.  **Migration:** Prefer MSCC-local helpers in `migration_helpers.R`
    over depending on package `MSdev`.
8.  **Demo molecule:** Glycine SMILES `NCC(O)=O` is the default in
    several helpers.

## External stacks

| Layer           | Packages / tools                              |
|-----------------|-----------------------------------------------|
| Core            | tidyverse, data.table, magrittr, stringr      |
| Structure       | ChemmineR (SDF), igraph, S4Vectors            |
| Spectra         | Spectra                                       |
| Optional chem   | enviPat, lc8, CHNOSZ, accucor, nnls           |
| Optional Python | reticulate → RDKit, RXNMapper                 |
| Optional UI     | shiny, plotly, visNetwork                     |
| Optional MS     | MSdev (Suggests only)                         |
| CFM-ID          | Docker container (see `CFM_*` in `dev_CFM.R`) |

## Docs & tests

- Vignettes: `vignettes/MSCC.Rmd` (stub), `vignettes/CFM_shiny.qmd` (CFM
  viewer usage).
- Informal tests/benchmarks: `Test/TS.R` (chemform_mz scaling).
- Man pages under `man/` are roxygen-generated from `R/`.

## High-value entry points

| Task | Start here |
|----|----|
| Formula → m/z | [`chemform_mz()`](https://drruili.github.io/MSCC/reference/chemform_mz.md), [`chemform_parse()`](https://drruili.github.io/MSCC/reference/chemform_parse.md) |
| Formula + adduct | [`chemform_adduct()`](https://drruili.github.io/MSCC/reference/chemform_adduct.md), [`chemform_adduct_check()`](https://drruili.github.io/MSCC/reference/chemform_adduct_check.md) |
| Isotope pattern | [`chemform_isotopes_pattern_enviPat()`](https://drruili.github.io/MSCC/reference/chemform_isotopes_pattern_enviPat.md) |
| Molecule graph | `get_Molecule_igraph_from_smiles()` |
| CFM predict | [`get_CFM_data_from_smiles()`](https://drruili.github.io/MSCC/reference/get_CFM_data_from_smiles.md) / [`CFM_predict()`](https://drruili.github.io/MSCC/reference/CFM.md) |
| CFM UI | [`shiny_vis_cfm()`](https://drruili.github.io/MSCC/reference/shiny_vis_cfm.md) |
| RDKit | [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md), [`rdkit_mol_formula()`](https://drruili.github.io/MSCC/reference/rdkit_mol_formula.md) |
| Natural correction | [`accucor_natural_correction()`](https://drruili.github.io/MSCC/reference/accucor_natural_correction.md) |
| Atom mapping | [`get_atom_map()`](https://drruili.github.io/MSCC/reference/get_atom_map.md), [`get_Molecule_atom_transfer_by_atom_map()`](https://drruili.github.io/MSCC/reference/Molecule_atom_transfer.md), `RXNMapper_map()` |

## Agent guardrails

- Treat this file as **read-only** unless the user explicitly asks to
  update it.
- Prefer editing exported APIs and their roxygen docs together;
  regenerate NAMESPACE via roxygen, do not hand-edit.
- Keep formula/isotope parsing behavior stable; many downstream MS
  workflows depend on string formats.
- `dev_*` files mix production exports and experimental code — check
  `@export` / NAMESPACE before assuming public API.
