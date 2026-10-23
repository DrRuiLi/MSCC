# Molecule graph from a SMILES string

Builds an igraph from the canonical RDKit molecule. Hydrogen atoms are
not included.

## Usage

``` r
get_sdf_igraph(sdf, addH = FALSE)
```

## Arguments

- sdf:

  A SMILES string, or several SMILES strings.

- addH:

  Ignored. Hydrogens are not part of the graph.

## Value

An igraph object, or a list of them when `sdf` has several strings.
