# Plot an RDKit molecule with SMARTS highlights

Same Cairo drawing as
[`rdkit_plot_mol()`](https://drruili.github.io/MSCC/reference/rdkit_plot_mol.md),
but fills atoms and bonds that match `highlight_smarts`. Named roles
`structure` / `co_elute` / `product` use gold, blue, and green.

## Usage

``` r
rdkit_plot_mol_highlight(
  mol = rdkit_mol_from_smiles(),
  highlight_smarts = NULL,
  width = 800L,
  height = 800L,
  bondLineWidth = 7,
  fixedFontSize = 80L,
  minFontSize = 72L,
  maxFontSize = 88L,
  padding = 0.12,
  additionalAtomLabelPadding = 0.04,
  multipleBondOffset = 0.16,
  scaleBondWidth = FALSE,
  centreMoleculesBeforeDrawing = TRUE,
  background = c(1, 1, 1)
)
```

## Arguments

- mol:

  RDKit molecule. Defaults to glycine via
  [`rdkit_mol_from_smiles()`](https://drruili.github.io/MSCC/reference/rdkit_mol_from_smiles.md).

- highlight_smarts:

  SMARTS string or character vector. Matching atoms and bonds are
  highlighted. Named values `structure` / `co_elute` / `product` select
  the palette. Unnamed values use gold.

- width, height:

  Integer pixel size of the Cairo canvas. Default `800`.

- bondLineWidth:

  Bond pen width. Default `7`.

- fixedFontSize:

  Atom-label size in pixels. Use `-1` to let RDKit scale with the
  molecule. Default `80`.

- minFontSize, maxFontSize:

  Font-size bounds when `fixedFontSize` is unused. Defaults `72` and
  `88`.

- padding:

  Fraction of the canvas left as margin. Default `0.12`.

- additionalAtomLabelPadding:

  Extra gap around heteroatom labels. Default `0.04`.

- multipleBondOffset:

  Offset between lines of a double/triple bond. Default `0.16`.

- scaleBondWidth:

  If `TRUE`, bond width scales with molecule size. Default `FALSE`.

- centreMoleculesBeforeDrawing:

  Centre the drawing in the canvas. Default `TRUE`.

- background:

  RGB vector in `[0, 1]`. Length 3 or 4 (alpha). Default white.

## Value

A ggplot object with the molecule raster.

## Examples

``` r
if (FALSE) { # \dontrun{
mol <- rdkit_mol_from_smiles("C(=O)(N)N")
rdkit_plot_mol_highlight(mol, c(structure = "[NX3H2][CX3]=O"))
} # }
```
