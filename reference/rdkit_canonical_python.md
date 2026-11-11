# RDKit interpreter that matches ICEBERG.

RDKit interpreter that matches ICEBERG.

## Usage

``` r
rdkit_canonical_python(python = NULL)
```

## Arguments

- python:

  Optional path. When omitted, uses `MSCC_RDKIT_PYTHON` if set,
  otherwise the ICEBERG CPU interpreter.
  [`ensure_RDKit_python()`](https://drruili.github.io/MSCC/reference/ensure_RDKit_python.md)
  selects the same executable for reticulate.

## Value

The python executable.
