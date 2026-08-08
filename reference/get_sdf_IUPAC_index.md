# Map ChemmineR SDF atom IDs to openclatura IUPAC locants

Converts an SDF to an in-memory molblock, loads it with RDKit
(preserving atomblock order), runs openclatura numbering, and returns a
table aligning ChemmineR canonical atom IDs with parent-chain IUPAC
locants.

## Usage

``` r
get_sdf_IUPAC_index(
  sdf,
  condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
  python = getOption("MSCC.rdkit.python", NULL)
)
```

## Arguments

- sdf:

  A ChemmineR `SDF` or `SDFset`.

- condaenv, python:

  Passed to
  [`ensure_RDKit_python()`](https://drruili.github.io/MSCC/reference/ensure_RDKit_python.md)
  before import.

## Value

For an `SDF`: a data.frame with columns `Atom_id`, `element`,
`rdkit_idx`, `locant`, `IUPAC_id`, and attribute `iupac_name`. For an
`SDFset`: a named list of such tables.

## Details

Requires a Python env with `rdkit` and `openclatura` (see
[`ensure_RDKit_python()`](https://drruili.github.io/MSCC/reference/ensure_RDKit_python.md)).

## Limitations

Locants come only from openclatura's **parent-skeleton NUMBERING** step
(`atom_to_locant`). Atoms that are not on that parent (substituents,
heteroatoms off the chain, most ring atoms when another fragment is
chosen as parent) stay `NA`. This is openclatura behavior, not an
SDF/RDKit index mismatch.

Example: for a steroid acetate named as *...-yl acetate*, openclatura
treats the acetate as the parent, so only the two acetate carbons
receive locants (`C1`/`C2`); the tetracyclic carbons remain `NA`. Stereo
locants embedded in the name string (e.g. `1S,2R,11S,...`) are not
exported as a per-atom map. Biological numbering (e.g. steroid C1–C17)
is a different nomenclature and is not provided here.

## Examples

``` r
if (FALSE) { # \dontrun{
sdf <- get_smiles_sdf("NCC(O)=O")[[1]]
idx <- get_sdf_IUPAC_index(sdf)
# C_2 -> C2 (alpha), C_3 -> C1 (carboxyl); N/O often NA
attr(idx, "iupac_name")
} # }
```
