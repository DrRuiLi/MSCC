# Test whether a mol matches a SMARTS pattern

Test whether a mol matches a SMARTS pattern

## Usage

``` r
rdkit_has_substruct(
  mol = rdkit_mol_from_smiles(),
  smarts = "[NH2]",
  Chem = get_RDKit_Chem()
)
```

## Arguments

- mol:

  RDKit molecule. Defaults to glycine via
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).

- smarts:

  Character SMARTS pattern. Default `"[NH2]"` (primary amine).

- Chem:

  Optional `rdkit.Chem` module from
  [`get_RDKit_Chem()`](https://drruili.github.io/MSCC/reference/get_RDKit_Chem.md).

## Value

Logical scalar.

## Examples

``` r
if (FALSE) { # \dontrun{
rdkit_has_substruct()
rdkit_has_substruct(rdkit_mol_from_smiles("CCO"), "[OH]")
} # }
```
