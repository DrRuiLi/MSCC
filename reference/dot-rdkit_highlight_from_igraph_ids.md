# Highlight atoms by Molecule_igraph vertex ids / map probabilities.

Highlight atoms by Molecule_igraph vertex ids / map probabilities.

## Usage

``` r
.rdkit_highlight_from_igraph_ids(mol, mol.ig, highlight_ids)
```

## Arguments

- mol:

  An RDKit molecule matching `mol.ig` atom order.

- mol.ig:

  A `Molecule_igraph`.

- highlight_ids:

  Character atom ids, or named numeric probabilities.

## Value

Highlight list for `.rdkit_draw_mol_ggplot()`.
