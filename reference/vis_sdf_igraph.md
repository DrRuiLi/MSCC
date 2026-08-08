# Visualize SDF igraph

Creates an interactive visualization of a molecular igraph using
visNetwork.

## Usage

``` r
vis_sdf_igraph(sdf.igraph, show_id = F, ...)
```

## Arguments

- sdf.igraph:

  An igraph object representing the molecule.

- show_id:

  Logical, whether to show atom IDs as labels (default FALSE).

- ...:

  Additional arguments passed to visNetwork.

## Value

A visNetwork html widget.
