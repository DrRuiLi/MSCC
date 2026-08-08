# RDKit bridge via reticulate

Helpers for importing RDKit submodules and building molecules from
SMILES. Import individual RDKit submodules by name (e.g.
[`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md),
or `reticulate::import("rdkit.Chem.Descriptors")`) — do not use a single
catch-all `get_RDKit()` for the whole package.

## Details

Requires a Python environment with `rdkit` configured for reticulate. By
default helpers select the conda env named `env_for_r` (override with
option `MSCC.rdkit.condaenv` or argument `condaenv`). Call this *before*
any other reticulate Python init in the session.
