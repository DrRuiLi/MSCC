# RDKit bridge via reticulate

Helpers for importing RDKit submodules and building molecules from
SMILES. Import individual RDKit submodules by name (e.g.
[`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md),
or `reticulate::import("rdkit.Chem.Descriptors")`) — do not use a single
catch-all `get_RDKit()` for the whole package.

## Details

Requires a Python environment with `rdkit` configured for reticulate.
Call
[`rdkit_reticulate_python`](https://drruili.github.io/MSCC/reference/rdkit_reticulate_python.md)
*before* any other reticulate Python init in the session. Default
interpreter: `MSCC_RDKIT_RETICULATE_PYTHON`, then
`MSCC_RDKIT_CLI_PYTHON`, then the ICEBERG CPU path (see
[`rdkit_cli_python`](https://drruili.github.io/MSCC/reference/rdkit_cli_python.md)).
CLI / subprocess helpers use
[`rdkit_cli_python`](https://drruili.github.io/MSCC/reference/rdkit_cli_python.md)
separately.
