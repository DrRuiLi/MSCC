# Select the Python env used for RDKit (before first import)

[`reticulate::import()`](https://rstudio.github.io/reticulate/reference/import.html)
cannot take an env path; the interpreter must be chosen first. On
Windows, also adds the conda env `Library/bin` directory to the DLL
search path so `rdMolDraw2D` / Cairo drawing can load. That add is *once
per Python process*: each `os.add_dll_directory()` cookie is retained,
because discarding them leaks Win32 search entries until `WinError 206`
(the path itself is not too long).

## Usage

``` r
ensure_RDKit_python(condaenv = NULL, python = NULL)
```

## Arguments

- condaenv:

  Optional conda env name or path. Used only when `python` is not set.
  The default is `NULL`, which selects
  [`rdkit_canonical_python`](https://drruili.github.io/MSCC/reference/rdkit_canonical_python.md).

- python:

  Optional path to a `python` executable. If set, takes precedence over
  `condaenv` and `MSCC_RDKIT_PYTHON`.

## Value

`invisible(TRUE)`.
