# Python executable for RDKit CLI / subprocess calls

Resolves the interpreter used by `system2` for
`inst/rdkit/canonical_mol.py` (canonical SMILES, molecule parts, MCS,
rxn ids). This is not the reticulate binder; see
[`rdkit_reticulate_python`](https://drruili.github.io/MSCC/reference/rdkit_reticulate_python.md).

## Usage

``` r
rdkit_cli_python(python = NULL)
```

## Arguments

- python:

  Optional path. When omitted, uses `MSCC_RDKIT_CLI_PYTHON` if set,
  otherwise the ICEBERG CPU interpreter.

## Value

The python executable path.
