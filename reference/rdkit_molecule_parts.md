# Atoms, bonds, formula, and a temporary molblock for one SMILES.

The molblock is only for callers that must hand a file to another
program. It is not stored on a molecule object.

## Usage

``` r
rdkit_molecule_parts(smiles)
```

## Arguments

- smiles:

  One SMILES string.

## Value

A list with smiles, formula, molblock, atoms, and bonds.
