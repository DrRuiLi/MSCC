# Map atoms between parent and product molecules

Performs atom mapping between two molecular structures using maximum
common substructure (MCS) analysis. Computes probability mappings of
atoms from parent to product, accounting for ring differences and bond
type similarities.

## Usage

``` r
get_atom_map(
  sdf.parent,
  sdf.product,
  ig.parent = get_sdf_igraph(sdf.parent),
  ig.product = get_sdf_igraph(sdf.product),
  iso_ele = "[13]C",
  return.type = c("most_prob", "prob_matrix")
)
```

## Arguments

- sdf.parent:

  An SDF object representing the parent molecule.

- sdf.product:

  An SDF object representing the product molecule.

- ig.parent:

  An igraph object representing the parent molecule's molecular graph.
  If NULL, it is computed from `sdf.parent`.

- ig.product:

  An igraph object representing the product molecule's molecular graph.
  If NULL, it is computed from `sdf.product`.

- iso_ele:

  Character string specifying the isotope element to consider for
  mapping (default: `"[13]C"`).

- return.type:

  Character string indicating the type of mapping to return. Either
  "most_prob" (default) returns a vector of most likely atom mappings,
  or "prob_matrix" returns a probability matrix.

## Value

For `return.type = "most_prob"`: a vector where names are atoms of the
product and values are probabilities of mapping to atoms of the parent.
For `return.type = "prob_matrix"`: a matrix with rows as atoms of the
parent and columns as atoms of the product, containing probabilities.
Both return values have an attribute "bond.score" indicating the bond
similarity score.
