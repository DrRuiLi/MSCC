# Plot an RDKit molecule (in-memory Draw)

Renders `mol` with RDKit `Draw.MolToImage` entirely in memory (no temp
file) and returns a ggplot that prints like a normal plot.

## Usage

``` r
rdkit_plot_mol(mol = rdkit_mol_from_smiles(), width = 500, height = 500L)
```

## Arguments

- mol:

  RDKit molecule. Defaults to glycine via
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).

- width, height:

  Integer pixel size passed to RDKit Draw.

## Value

A ggplot object with the molecule raster.

## Examples

``` r
if (FALSE) { # \dontrun{
rdkit_plot_mol()
rdkit_plot_mol(rdkit_mol_from_smiles("CCO"))
} # }
```
