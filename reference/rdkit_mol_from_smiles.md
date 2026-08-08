# Build an RDKit molecule from SMILES

Build an RDKit molecule from SMILES

## Usage

``` r
rdkit_mol_from_smiles(smiles = "NCC(O)=O", Chem = get_RDKit_Chem())
```

## Arguments

- smiles:

  Character scalar SMILES. Default `"NCC(O)=O"` (glycine), matching MSCC
  helpers such as `get_Molecule_igraph_from_smiles()`.

- Chem:

  Optional `rdkit.Chem` module from
  [`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md).

## Value

An RDKit mol object, or `NULL` if invalid.

## Examples

``` r
if (FALSE) { # \dontrun{
mol <- rdkit_mol_from_smiles()
mol <- rdkit_mol_from_smiles("CCO")
} # }
```
