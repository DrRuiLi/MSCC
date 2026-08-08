# Predict Mass Spectra using CFM-ID

Predicts mass spectra fragments for a given molecule using the CFM-ID
algorithm via Docker. This function sends a SMILES or InChI string to
the CFM-ID Docker container and returns predicted fragment spectra at
multiple collision energies.

Annotates experimental mass spectra with predicted fragments using the
CFM-ID algorithm via Docker. This function matches experimental spectra
to predicted fragments based on mass tolerance.

Deprecated. Use `CFM_annotate_by_predict` instead. This function
generates CFM_data by generating fragment graphs with CFM_fraggen, not
relying on Spectra input. It creates a synthetic peak assignment based
on the generated fragments.

Combines `CFM_predict` and `CFM_annotate` in a single workflow. First
predicts mass spectra using CFM-ID, then annotates the predicted spectra
with fragment assignments. This is useful for obtaining fully annotated
spectral data for a given molecule.

Optional disk cache (`check_cache`): stores/loads
`<id>_positive_CFM_data.rds` / `<id>_negative_CFM_data.rds` under
`cache_dir`.

Reads and parses the output from CFM-ID prediction results into a
structured CFM_data object. This function processes the text output file
from CFM-ID and extracts peak assignments, fragment definitions, and
other spectral data.

## Usage

``` r
CFM_predict(
  smiles_or_inchi_or_file = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
  prob_thresh = 0.001,
  param_adduct = "[M+H]+",
  annotate_fragments = 1,
  output_file_or_dir = NULL,
  apply_postproc = 0,
  suppress_exceptions = 1
)

CFM_annotate(
  smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
  spectrum_file = NULL,
  id = "AN_ID",
  ppm_mass_tol = 5,
  abs_mass_tol = 0,
  param_adduct = "[M+H]+",
  output_file = NULL,
  ...
)

CFM_annotate_by_fraggen(
  smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
  spectrum_file = NULL,
  max_depth = 1,
  id = "AN_ID",
  ppm_mass_tol = 5,
  abs_mass_tol = 0,
  param_adduct = "[M+H]+",
  output_file = NULL,
  ...
)

CFM_annotate_by_predict(
  smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
  id = NULL,
  ppm_mass_tol = 5,
  abs_mass_tol = 0,
  param_adduct = "[M+H]+",
  output_file = NULL,
  check_cache = FALSE,
  cache_dir = tempdir(),
  force = FALSE,
  ...
)

CFM_fraggen(
  smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
  max_depth = 2,
  param_adduct = "[M+H]+",
  output_file_or_dir = NULL
)

read_CFM_predict_result(
  result_path = "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_predict_result.txt",
  polarity = NA_real_
)
```

## Arguments

- smiles_or_inchi_or_file:

  SMILES string, InChI string, or path to a file containing molecular
  structure

- prob_thresh:

  Probability threshold for fragment prediction (default: 0.001)

- param_adduct:

  Adduct type, e.g., "\[M+H\]+" or "\[M-H\]-" (default: "\[M+H\]+").
  Also accepts `0`/`1` (0 = "\[M-H\]-", 1 = "\[M+H\]+") and `"-"`/`"+"`.

- annotate_fragments:

  Logical, whether to annotate fragments (default: 1)

- output_file_or_dir:

  Path to save results, or NULL to return results in memory (default:
  NULL)

- apply_postproc:

  Logical, whether to apply post-processing (default: 0)

- suppress_exceptions:

  Logical, whether to suppress exceptions (default: 1)

- smiles_or_inchi:

  SMILES string or InChI string of the molecule

- spectrum_file:

  Spectrum file path (currently unused, default: NULL)

- id:

  Compound identifier. If `NULL` or empty, a random id is generated each
  call (cache will not be reusable across runs unless `id` is supplied).

- ppm_mass_tol:

  Mass tolerance in ppm for annotation matching (default: 5.0)

- abs_mass_tol:

  Absolute mass tolerance in m/z units (default: 0)

- output_file:

  Path to save results, or NULL to return results in memory (default:
  NULL)

- ...:

  Additional arguments passed to underlying functions

- max_depth:

  Maximum depth for fragment generation (default: 1)

- check_cache:

  Logical; if TRUE, reuse/save RDS cache under `cache_dir`

- cache_dir:

  Cache directory (used when `check_cache` is TRUE). Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- force:

  Logical; if TRUE, ignore existing cache and recompute

- result_path:

  Path to the CFM-ID prediction result file (default:
  "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_predict_result.txt")

- polarity:

  Numeric polarity value (0 for negative, 1 for positive, default: NA)

## Value

A CFM_data object containing predicted spectra data, or path to output
file if output_file_or_dir is specified

A CFM_data object containing annotated spectra data, or path to output
file if output_file is specified

A CFM_data object containing fragment definitions and peak assignments

A CFM_data object containing parsed peak assignments and fragment
definitions

## Details

Predict Mass Spectra using CFM-ID

Annotate Mass Spectra using CFM-ID

Annotate Fragments using CFM Fraggen (Deprecated)

Read CFM-ID Prediction Results

## Functions

- `CFM_predict()`: predict

- `CFM_annotate()`: annotate

- `CFM_annotate_by_fraggen()`: fraggen and annotate

- `CFM_annotate_by_predict()`: predict and annotate

- `CFM_fraggen()`: fraggen

- `read_CFM_predict_result()`: read prediction results

## Note

Filters fragments with mz \< 30 to avoid errors when converting SMILES
to SDF. This function is deprecated and will be removed in future
versions.
