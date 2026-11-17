# Bind reticulate to the RDKit Python interpreter (before first import)

[`reticulate::import()`](https://rstudio.github.io/reticulate/reference/import.html)
cannot take an env path; the interpreter must be chosen first. On
Windows, also adds the conda env `Library/bin` directory to the DLL
search path so `rdMolDraw2D` / Cairo drawing can load. That add is *once
per Python process*: each `os.add_dll_directory()` cookie is retained,
because discarding them leaks Win32 search entries until `WinError 206`
(the path itself is not too long).

## Usage

``` r
rdkit_reticulate_python(condaenv = NULL, python = NULL)
```

## Arguments

- condaenv:

  Optional conda env name or path. Used only when `python` is not set.

- python:

  Optional path to a `python` executable. If set, takes precedence over
  `condaenv` and the env vars above.

## Value

`invisible(TRUE)`.

## Details

Resolution when `python` and `condaenv` are both unset:
`MSCC_RDKIT_RETICULATE_PYTHON`, then `MSCC_RDKIT_CLI_PYTHON`, then the
ICEBERG hardcoded path via
[`rdkit_cli_python`](https://drruili.github.io/MSCC/reference/rdkit_cli_python.md).
