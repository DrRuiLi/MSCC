# Plot a SMILES string with ggplot geometry

Like
[`ggplot_Molecule_igraph()`](https://drruili.github.io/MSCC/reference/ggplot_Molecule_igraph.md),
but takes SMILES directly. Useful for CFM fragment SMILES that fail
strict sanitization / canonicalization (for example `C=CC#CC#[CH2+]`);
those are drawn with `sanitize = FALSE`.

## Usage

``` r
ggplot_smiles(
  smiles,
  atom_label = c("element", "id", "hetero"),
  sanitize = TRUE,
  highlight_ids = NULL,
  label_size = 3.2,
  bond_linewidth = 0.7,
  highlight_size = 8
)
```

## Arguments

- smiles:

  Character scalar SMILES.

- atom_label:

  `"element"` labels every atom with its element (and formal charge when
  non-zero), `"id"` uses `C_1`-style ids, `"hetero"` hides neutral
  carbons.

- sanitize:

  Passed to
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).
  If `TRUE` and parsing fails, automatically retries with
  `sanitize = FALSE`.

- highlight_ids:

  Optional atom ids / named probabilities.

- label_size, bond_linewidth, highlight_size:

  Drawing sizes.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot_smiles("CCO", atom_label = "element")
ggplot_smiles("C=CC#CC#[CH2+]", atom_label = "element")
} # }
```
