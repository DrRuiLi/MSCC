# Plot a Molecule_igraph with ggplot geometry

Reads 2D coordinates, Kekule bond orders, ring membership, and wedge
directions from RDKit after `PrepareMolForDrawing()`. Bonds are
`geom_segment` layers and atom ids are `geom_text`, so a PDF device
keeps them as vectors. Element colours follow RDKit's default palette. A
ring double bond is offset toward the ring centre; a non-ring double
bond is a symmetric pair of lines.

## Usage

``` r
ggplot_Molecule_igraph(mol.ig, show_id = TRUE)
```

## Arguments

- mol.ig:

  A `Molecule_igraph`.

- show_id:

  If `TRUE`, label every atom with its vertex id. If `FALSE`, label
  heteroatoms with the element symbol. Default `TRUE`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
mig <- get_Molecule_igraph_from_smiles("NCC(O)=O")
ggplot_Molecule_igraph(mig)
} # }
```
