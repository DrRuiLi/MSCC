# AI_CONTEXT — MSCC

> Architectural map for AI assistants. Read this before exploring the
> codebase. Do **not** auto-update this file; propose changes and wait
> for explicit approval.

## Package identity

| Field | Value |
|----|----|
| Name | **MSCC** (Mass Spectrometry Chemical Calculation) |
| Version | 1.3.1 |
| Author | Rui Li (`rli@sinh.ac.cn`) |
| Repo | <https://github.com/drruili/MSCC> |
| Type | R package (S4 + tidyverse + igraph + RDKit; optional reticulate / Spectra) |

**Purpose:** Chemical-formula arithmetic, isotope/adduct handling for
MS, molecular graph representation, CFM-ID spectrum
prediction/annotation, and related metabolomics helpers (atom mapping,
natural-abundance correction, interactive viewers).

## Mental model

    SMILES / formula / SDF
            │
            ├─► chemform_*  ──► m/z, adducts, isotope patterns
            │
            ├─► Molecule_igraph (RDKit CLI graph + igraph + isotopomers)
            │         │
            │         ├─► atom maps / MCS (CLI); RXNMapper deprecated
            │         └─► rdkit_plot_* / visNetwork
            │
            └─► CFM_data (Docker CFM-ID)
                      │
                      ├─► Spectra / peak lists
                      └─► shiny_vis_cfm / plotly_CFM_spectra

## Core S4 classes

### `Molecule_igraph` (`R/Molecule_igraph-class.R`)

Slots: `molecule_info`, `sdf` (legacy ChemmineR `SDF` where still used),
`igraph`, `isotopomer` (data.frame).

- Build via `get_Molecule_igraph_from_smiles()` /
  `get_Molecule_igraph_from_sdf()` (graphs from RDKit CLI
  [`rdkit_molecule_parts()`](https://drruili.github.io/MSCC/reference/rdkit_molecule_parts.md)).
- Graph attrs via generics `vdata` / `edata` (defined in
  `R/00_graph_primitives.R`).
- Isotopomer helpers:
  [`Molecule_igraph_add_isotopomer()`](https://drruili.github.io/MSCC/reference/Molecule_igraph_add_isotopomer.md);
  atom transfer via MCS (`rdkit_mcs_pairs`). RXNMapper path is
  deprecated.

### `CFM_data` (`R/dev_CFM.R`)

Slots:

| Slot | Role |
|----|----|
| `peak_assignment` | Peaks per CE (`energy0/1/2` → CE 10/20/40), optional `fragment_id` / `fragment_score` |
| `fragment_define` | Fragment SMILES, `fragment_mz`, polarity |
| `fragment_transition` | Fragment graph edges (annotate path) |

Typical entry: `get_CFM_data_from_smiles()` → `CFM_predict` /
`CFM_annotate*` (Docker). Viewer: `shiny_vis_cfm()`.

## Module map (`R/`)

| File | Responsibility |
|----|----|
| `Chemform-function.R` | Parse/format formulas; `chemform_parse`, `chemform_get_ele`, `chemform_formate`, `chemform_mz`, `chemform_calc` |
| `Chemform_calc.R` | Formula sum/multiply (`chemform_sum`, `chemform_multi`) |
| `Adduct-function.R` | Adduct table logic; `chemform_adduct`, `chemform_adduct_check` |
| `isotope.R` | Isotope mass diffs, formula labeling (`chemform_isotope_label`) |
| `dev_envipat_lc8.R` | enviPat isotope patterns (`chemform_isotopes_pattern_enviPat`) |
| `migration_helpers.R` | Self-contained helpers (no MSdev); elem table, xcms isotope matching, polarity |
| `accucor_wrapper.R` | Natural-isotope correction (`accucor_natural_correction`) |
| `00_graph_primitives.R` | `vdata`/`edata`/`atom`/`get_element` generics; deprecated `get_RXNMapper` |
| `Molecule_igraph-class.R` | S4 class + methods |
| `Molecule_igraph-functions.R` | Build/filter/vis isotopomers; atom-transfer matrices |
| `Molecule_vis.R` | SDF↔︎igraph, `vis_sdf_igraph` |
| `dev_Chemmine.R` | Legacy SMILES↔︎SDF / vis helpers; prefer RDKit CLI entries |
| `dev_CFM.R` | CFM-ID Docker API, `CFM_data` I/O, Spectra bridges |
| `shiny_vis_cfm.R` | Shiny + plotly CFM viewer |
| `rdkit_canonical.R` | CLI Python path + `system2` → `inst/rdkit/canonical_mol.py` |
| `dev_rdkit.R` | reticulate RDKit (`rdkit_reticulate_python`, `get_RDKit_Chem`, plot/SMARTS) |
| `dev_RXNmapper.R` | Deprecated RXNMapper wrappers (`RXNMapper_map`, …) |
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
5.  **RDKit Python (two paths):**
    - **CLI / subprocess:**
      [`rdkit_cli_python()`](https://drruili.github.io/MSCC/reference/rdkit_cli_python.md)
      → `system2` + `inst/rdkit/canonical_mol.py`. Env:
      `MSCC_RDKIT_CLI_PYTHON`, else ICEBERG
      `/mnt/nas/conda/envs/lirui_iceberg_cpu/bin/python`.
    - **reticulate:**
      [`rdkit_reticulate_python()`](https://drruili.github.io/MSCC/reference/rdkit_reticulate_python.md)
      before first import. Env: `MSCC_RDKIT_RETICULATE_PYTHON` →
      `MSCC_RDKIT_CLI_PYTHON` → same ICEBERG path. Optional `condaenv=`
      / `python=` overrides.
    - Import RDKit submodules explicitly
      ([`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md));
      do not invent a catch-all `get_RDKit()`. Old names
      `rdkit_canonical_python` / `ensure_RDKit_python` /
      `MSCC_RDKIT_PYTHON` are gone.
6.  **RXNMapper:** `get_RXNMapper` / `RXNMapper_map` /
    `RXNMapper_mapped_rxn_parse` are **deprecated** (still exported).
    Prefer
    [`rdkit_mcs_pairs()`](https://drruili.github.io/MSCC/reference/rdkit_mcs_pairs.md)
    /
    [`rdkit_rxn_ids()`](https://drruili.github.io/MSCC/reference/rdkit_rxn_ids.md)
    for now.
7.  **CFM runtime:** Prediction/annotation needs a working CFM-ID Docker
    setup; local cache dirs are supported.
8.  **No MSdev dependency:** Do **not** rely on package `MSdev` (no
    `Imports`/`Suggests`/`Remotes`, no `MSdev::` calls). MSdev already
    depends on MSCC; adding the reverse creates a circular install loop.
    Use MSCC-local helpers in `migration_helpers.R` and CFM plotting via
    `plotly_CFM_spectra()` / `shiny_vis_cfm()`.
9.  **Demo molecule:** Glycine SMILES `NCC(O)=O` is the default in
    several helpers.

## External stacks

| Layer | Packages / tools |
|----|----|
| Core | tidyverse, data.table, magrittr, stringr |
| Structure | ChemmineR (SDF), igraph, S4Vectors |
| Spectra | Spectra |
| Optional chem | enviPat, CHNOSZ, accucor, nnls |
| Optional Python | CLI RDKit (`rdkit_cli_python`); reticulate RDKit (`rdkit_reticulate_python`); RXNMapper deprecated |
| Optional UI | shiny, plotly, visNetwork |
| Downstream (depends on MSCC) | MSdev — one-way only; MSCC must not depend back |
| CFM-ID | Docker container (see `CFM_*` in `dev_CFM.R`) |

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
| CFM predict | `get_CFM_data_from_smiles()` / `CFM_predict()` |
| CFM UI | `shiny_vis_cfm()` |
| RDKit CLI | [`rdkit_cli_python()`](https://drruili.github.io/MSCC/reference/rdkit_cli_python.md), [`rdkit_molecule_parts()`](https://drruili.github.io/MSCC/reference/rdkit_molecule_parts.md), [`rdkit_canonical_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_canonical_smiles.md) |
| RDKit reticulate | [`rdkit_reticulate_python()`](https://drruili.github.io/MSCC/reference/rdkit_reticulate_python.md), [`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md), [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md), [`rdkit_plot_mol()`](https://drruili.github.io/MSCC/reference/rdkit_plot_mol.md) |
| Natural correction | [`accucor_natural_correction()`](https://drruili.github.io/MSCC/reference/accucor_natural_correction.md) |
| Atom mapping | [`rdkit_mcs_pairs()`](https://drruili.github.io/MSCC/reference/rdkit_mcs_pairs.md), [`rdkit_rxn_ids()`](https://drruili.github.io/MSCC/reference/rdkit_rxn_ids.md) (RXNMapper APIs deprecated) |

## Agent guardrails

- Treat this file as **read-only** unless the user explicitly asks to
  update it.
- Prefer editing exported APIs and their roxygen docs together;
  regenerate NAMESPACE via roxygen, do not hand-edit.
- Keep formula/isotope parsing behavior stable; many downstream MS
  workflows depend on string formats.
- `dev_*` files mix production exports and experimental code — check
  `@export` / NAMESPACE before assuming public API.
- **Never add `MSdev` to DESCRIPTION or call `MSdev::` from MSCC** —
  that creates a circular dependency (MSdev → MSCC already).
