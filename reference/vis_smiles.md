# Visualize a molecule from SMILES string

Creates an interactive visualization of a molecular structure from a
SMILES string. Converts the SMILES to an SDF object, then to an igraph
representation, and generates an HTML widget using visNetwork.
Optionally displays the molecular formula.

## Usage

``` r
vis_smiles(smiles, show.formula = T, show_id = T, highlight = NULL)
```

## Arguments

- smiles:

  A single SMILES string representing the molecule to visualize.

- show.formula:

  Logical indicating whether to display the molecular formula in the
  plot (default: TRUE).

- show_id:

  Logical indicating whether to show atom IDs as labels (default: TRUE).
  If FALSE, atom symbols are shown.

- highlight:

  Optional character vector of atom IDs to highlight in the
  visualization.

## Value

An HTML widget object (visNetwork) that can be rendered in RStudio
viewer or browser.
