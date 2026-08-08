# Molecular formula from an RDKit mol

Molecular formula from an RDKit mol

## Usage

``` r
rdkit_mol_formula(mol = rdkit_mol_from_smiles())
```

## Arguments

- mol:

  RDKit molecule. Defaults to glycine via
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).

## Value

Character formula string, or `NA_character_`.

## Examples

``` r
if (FALSE) { # \dontrun{
rdkit_mol_formula()
rdkit_mol_formula(rdkit_mol_from_smiles("CCO"))
} # }
```
