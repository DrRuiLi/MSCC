# Add isotopomer to Molecule Igraph

Adds an isotopomer to a Molecule_igraph object, updating its isotopomer
data frame.

## Usage

``` r
Molecule_igraph_add_isotopomer(
  Molecule_igraph,
  isotopomer = NULL,
  iso_vec = NULL,
  abundance = NA,
  path = NA,
  FSIS = NA
)
```

## Arguments

- Molecule_igraph:

  A `Molecule_igraph` object.

- isotopomer:

  Isotopomer name. If NULL, a unique name is generated.

- iso_vec:

  Named vector mapping atom IDs to isotopic labels, e.g., c("C_1" =
  `"[13]C"`).

- abundance:

  Numeric abundance of the isotopomer.

- path:

  Path information (optional).

- FSIS:

  Fragment-specific isotopomer score (optional).

## Value

The updated `Molecule_igraph` object.
