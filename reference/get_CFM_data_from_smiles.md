# Get CFM Data from SMILES

Creates CFM_data object from SMILES by running CFM prediction and
annotation.

## Usage

``` r
get_CFM_data_from_smiles(
  smiles = "NCC(O)=O",
  compound_id = NULL,
  ppm = 5,
  adduct = "[M+H]+",
  check_cache = FALSE,
  cache_dir = tempdir(),
  force = FALSE,
  ...
)
```

## Arguments

- smiles:

  SMILES string of the molecule

- compound_id:

  Identifier for the compound. If `NULL` or empty, a random id is
  generated each call (see
  [`CFM_annotate_by_predict`](https://drruili.github.io/MSCC/reference/CFM.md)).

- ppm:

  Mass tolerance in ppm (default: 5)

- adduct:

  Adduct type (default: "\[M+H\]+"). Also accepts 0/1 (0 = "\[M-H\]-", 1
  = "\[M+H\]+") and '-'/'+' as shorthand.

- check_cache:

  Logical; forward to
  [`CFM_annotate_by_predict`](https://drruili.github.io/MSCC/reference/CFM.md)

- cache_dir:

  Cache directory; forward to
  [`CFM_annotate_by_predict`](https://drruili.github.io/MSCC/reference/CFM.md)
  (default [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- force:

  Logical; ignore cache if TRUE

- ...:

  Additional arguments

## Value

A CFM_data object containing fragment data

## Details

Get CFM Data from SMILES
