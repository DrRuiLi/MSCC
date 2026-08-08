# Convert SMILES strings to SDF format

Converts one or more SMILES strings to SDF (Structure Data Format)
objects using the ChemmineR package. Optionally canonicalizes the
structures. Uses a precomputed mapping table to replace known SMILES
with stored SDFs.

## Usage

``` r
get_smiles_sdf(smiles, smiles.id = names(smiles), canonicalize = T)
```

## Arguments

- smiles:

  Character vector of SMILES strings to convert.

- smiles.id:

  Optional character vector of IDs to assign to the resulting SDF
  objects. If NULL and `smiles` has names, those names are used;
  otherwise IDs are generated as "CMP001", etc.

- canonicalize:

  Logical indicating whether to canonicalize the SDF structures
  (default: TRUE).

## Value

An SDFset object (list of SDF objects) containing the molecular
structures.

## Details

get_smiles_sdf
