# Build highlight points for [`ggplot_Molecule_igraph()`](https://drruili.github.io/MSCC/reference/ggplot_Molecule_igraph.md).

Build highlight points for
[`ggplot_Molecule_igraph()`](https://drruili.github.io/MSCC/reference/ggplot_Molecule_igraph.md).

## Usage

``` r
.ggplot_highlight_atoms(atoms, highlight_ids, highlight_size = 8)
```

## Arguments

- atoms:

  Atom table from `.ggplot_molecule_drawing()`.

- highlight_ids:

  Character ids or named numeric probabilities.

- highlight_size:

  Point size in mm.

## Value

A data.frame with x, y, fill, size (possibly 0 rows).
