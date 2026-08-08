# Convert SDF to igraph Object

Converts a chemical structure from SDF format to an igraph object for
network visualization.

## Usage

``` r
get_sdf_igraph(sdf, addH = F)
```

## Arguments

- sdf:

  An SDF or SDFset object.

- addH:

  Logical, whether to include hydrogen atoms (default FALSE).

## Value

An igraph object representing the molecular structure.
