# Import the `rdkit.Chem` Python module

Returns only `rdkit.Chem`. For other RDKit submodules, import them
explicitly (e.g. `reticulate::import("rdkit.Chem.Descriptors")`).

## Usage

``` r
get_RDKit_Chem(
  condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
  python = getOption("MSCC.rdkit.python", NULL)
)
```

## Arguments

- condaenv, python:

  Passed to
  [`ensure_RDKit_python()`](https://drruili.github.io/MSCC/reference/ensure_RDKit_python.md)
  before import.

## Value

The imported `rdkit.Chem` module.
