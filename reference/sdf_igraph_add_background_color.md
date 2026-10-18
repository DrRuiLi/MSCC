# Color atom backgrounds on a structure graph.

Color atom backgrounds on a structure graph.

## Usage

``` r
sdf_igraph_add_background_color(
  sdf.igraph,
  value,
  color.ramp = colramp(breaks = c(0, Inf, 1), colors = c("#FFFFFF", "#F7844F",
    "#B20C26"))
)
```

## Arguments

- sdf.igraph:

  An igraph from a molecular structure.

- value:

  Named numeric weights. Names are atom ids.

- color.ramp:

  Color ramp used by the structure viewer.

## Value

The graph with `color.background` set on vertices.
