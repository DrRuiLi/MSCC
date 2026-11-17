# Plot a Molecule_igraph with RDKit atom ids

Draws `mol.ig@molecule_info$smiles` with
[`rdkit_plot_mol()`](https://drruili.github.io/MSCC/reference/rdkit_plot_mol.md).
When `show_id` is `TRUE`, each atom label is the graph vertex `id`
(`C_1`, `N_1`, ...). Those ids are the RDKit atom index plus one,
prefixed by the element, and they are the column names of `@isotopomer`.
The label is written as RDKit's `atomLabel` property before drawing,
which replaces the element symbol.

## Usage

``` r
rdkit_plot_Molecule_igraph(mol.ig, show_id = TRUE, highlight_ids = NULL, ...)
```

## Arguments

- mol.ig:

  A `Molecule_igraph`.

- show_id:

  If `TRUE`, label atoms with vertex ids. Default `TRUE`.

- highlight_ids:

  Optional atom ids to highlight, or a named numeric vector of mapping
  probabilities (names = vertex ids). Values `> 0` are highlighted;
  colour intensity follows the probability when names are supplied.
  Bonds between highlighted atoms are filled too.

- ...:

  Passed to
  [`rdkit_plot_mol()`](https://drruili.github.io/MSCC/reference/rdkit_plot_mol.md)
  / the Cairo drawer. With `show_id = TRUE`, `fixedFontSize` defaults to
  `28` and `padding` to `0.16` unless those arguments are supplied.

## Value

A ggplot object with the molecule raster.

## Details

Vertex row order must match the atom order of the canonical SMILES. The
function stops if the atom count or element symbols disagree.

## Examples

``` r
if (FALSE) { # \dontrun{
mig <- get_Molecule_igraph_from_smiles("NCC(O)=O")
rdkit_plot_Molecule_igraph(mig)
rdkit_plot_Molecule_igraph(mig, highlight_ids = c("C_2", "C_3"))
} # }
```
