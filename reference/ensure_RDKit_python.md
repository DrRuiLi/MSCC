# Select the Python env used for RDKit (before first import)

[`reticulate::import()`](https://rstudio.github.io/reticulate/reference/import.html)
cannot take an env path; the interpreter must be chosen first. On
Windows, also adds the conda env `Library/bin` directory to the DLL
search path so `rdMolDraw2D` / Cairo drawing can load.

## Usage

``` r
ensure_RDKit_python(
  condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
  python = getOption("MSCC.rdkit.python", NULL)
)
```

## Arguments

- condaenv:

  Conda env name or path. Default
  `getOption("MSCC.rdkit.condaenv", "env_for_r")`. Use `NULL` to skip
  auto-selection.

- python:

  Optional path to a `python` executable. If set, takes precedence over
  `condaenv`.

## Value

`invisible(TRUE)`.
