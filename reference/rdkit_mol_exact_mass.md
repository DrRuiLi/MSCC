# Monoisotopic exact mass from an RDKit mol

Monoisotopic exact mass from an RDKit mol

## Usage

``` r
rdkit_mol_exact_mass(mol = rdkit_mol_from_smiles())
```

## Arguments

- mol:

  RDKit molecule. Defaults to glycine via
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).

## Value

Numeric exact mass, or `NA_real_`.

## Examples

``` r
if (FALSE) { # \dontrun{
rdkit_mol_exact_mass()
rdkit_mol_exact_mass(rdkit_mol_from_smiles("CCO"))
} # }
```
