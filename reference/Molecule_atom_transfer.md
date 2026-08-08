# Auto map atom structure

Molecule atom transfer by atom map.

## Usage

``` r
get_Molecule_atom_transfer_by_atom_map(
  mol.ig.from,
  mol.ig.to,
  target_ele = "ANY"
)
```

## Arguments

- mol.ig.from:

  Molecule_igraph

- mol.ig.to:

  Molecule_igraph

## Value

`matirx`, column as atom of `to`, row as multiple map, value as atom of
`from`

## Functions

- `get_Molecule_atom_transfer_by_atom_map()`: return a matrix, column as
  atom of `to`, row as multiple map, value as atom of `from`
