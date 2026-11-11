# RDKit bridge via reticulate

Helpers for importing RDKit submodules and building molecules from
SMILES. Import individual RDKit submodules by name (e.g.
[`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md),
or `reticulate::import("rdkit.Chem.Descriptors")`) — do not use a single
catch-all `get_RDKit()` for the whole package.

## Details

Requires a Python environment with `rdkit` configured for reticulate.
Helpers use the same interpreter as
[`rdkit_canonical_python`](https://drruili.github.io/MSCC/reference/rdkit_canonical_python.md):
`MSCC_RDKIT_PYTHON` when set, otherwise
`/mnt/nas/conda/envs/lirui_iceberg_cpu/bin/python`. Call this *before*
any other reticulate Python init in the session.
