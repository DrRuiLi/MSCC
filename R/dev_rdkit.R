#' RDKit bridge via reticulate
#'
#' Helpers for importing RDKit submodules and building molecules from SMILES.
#' Import individual RDKit submodules by name (e.g. \code{get_RDKit_Chem()},
#' or \code{reticulate::import("rdkit.Chem.Descriptors")}) — do not use a
#' single catch-all \code{get_RDKit()} for the whole package.
#'
#' Requires a Python environment with \code{rdkit} configured for
#' \pkg{reticulate}. By default helpers select the conda env named
#' \code{env_for_r} (override with option \code{MSCC.rdkit.condaenv} or
#' argument \code{condaenv}). Call this \emph{before} any other reticulate
#' Python init in the session.
#'
#' @name dev_rdkit
NULL

#' Select the Python env used for RDKit (before first import)
#'
#' \code{reticulate::import()} cannot take an env path; the interpreter must
#' be chosen first. On Windows, also adds the conda env
#' \code{Library/bin} directory to the DLL search path so
#' \code{rdMolDraw2D} / Cairo drawing can load.
#'
#' @param condaenv Conda env name or path. Default
#'   \code{getOption("MSCC.rdkit.condaenv", "env_for_r")}. Use \code{NULL}
#'   to skip auto-selection.
#' @param python Optional path to a \code{python} executable. If set, takes
#'   precedence over \code{condaenv}.
#' @return \code{invisible(TRUE)}.
#' @export
ensure_RDKit_python <- function(condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
                                python = getOption("MSCC.rdkit.python", NULL)) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for RDKit. Install it and configure a Python env with rdkit.")
  }
  if (!isTRUE(reticulate::py_available(initialize = FALSE))) {
    if (!is.null(python) && nzchar(python)) {
      reticulate::use_python(python, required = TRUE)
    } else if (!is.null(condaenv) && nzchar(condaenv)) {
      reticulate::use_condaenv(condaenv, required = TRUE)
    }
  }
  .rdkit_add_windows_dll_dirs()
  invisible(TRUE)
}

# Windows Python >= 3.8 does not search PATH for extension DLL deps.
.rdkit_add_windows_dll_dirs <- function() {
  if (.Platform$OS.type != "windows") {
    return(invisible(FALSE))
  }
  if (!isTRUE(reticulate::py_available(initialize = TRUE))) {
    return(invisible(FALSE))
  }
  reticulate::py_run_string(
    paste(
      "import os",
      "from pathlib import Path",
      "import sys",
      "for p in (Path(sys.prefix) / 'Library' / 'bin', Path(sys.base_prefix) / 'Library' / 'bin'):",
      "    if p.is_dir():",
      "        os.add_dll_directory(str(p))",
      sep = "\n"
    )
  )
  invisible(TRUE)
}

#' Import the \code{rdkit.Chem} Python module
#'
#' Returns only \code{rdkit.Chem}. For other RDKit submodules, import them
#' explicitly (e.g. \code{reticulate::import("rdkit.Chem.Descriptors")}).
#'
#' @param condaenv,python Passed to [ensure_RDKit_python()] before import.
#' @return The imported \code{rdkit.Chem} module.
#' @export
get_RDKit_Chem <- function(condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
                           python = getOption("MSCC.rdkit.python", NULL)) {
  ensure_RDKit_python(condaenv = condaenv, python = python)
  reticulate::py_suppress_warnings(
    reticulate::import("rdkit.Chem", convert = TRUE)
  )
}

#' Build an RDKit molecule from SMILES
#'
#' @param smiles Character scalar SMILES. Default \code{"NCC(O)=O"} (glycine),
#'   matching MSCC helpers such as \code{get_Molecule_igraph_from_smiles()}.
#' @param Chem Optional \code{rdkit.Chem} module from [get_RDKit_Chem()].
#' @return An RDKit mol object, or \code{NULL} if invalid.
#'
#' @examples
#' \dontrun{
#' mol <- rdkit_mol_from_smiles()
#' mol <- rdkit_mol_from_smiles("CCO")
#' }
#'
#' @export
rdkit_mol_from_smiles <- function(smiles = "NCC(O)=O", Chem = get_RDKit_Chem()) {
  if (is.null(smiles) || length(smiles) != 1L || is.na(smiles) || !nzchar(smiles)) {
    return(NULL)
  }
  mol <- Chem$MolFromSmiles(as.character(smiles))
  if (reticulate::py_is_null_xptr(mol) || is.null(mol)) {
    return(NULL)
  }
  mol
}

#' Molecular formula from an RDKit mol
#'
#' @param mol RDKit molecule. Defaults to glycine via
#'   [rdkit_mol_from_smiles()].
#' @return Character formula string, or \code{NA_character_}.
#'
#' @examples
#' \dontrun{
#' rdkit_mol_formula()
#' rdkit_mol_formula(rdkit_mol_from_smiles("CCO"))
#' }
#'
#' @export
rdkit_mol_formula <- function(mol = rdkit_mol_from_smiles()) {
  if (is.null(mol)) {
    return(NA_character_)
  }
  desc <- reticulate::import("rdkit.Chem.rdMolDescriptors", convert = TRUE)
  as.character(desc$CalcMolFormula(mol))
}

#' Monoisotopic exact mass from an RDKit mol
#'
#' @param mol RDKit molecule. Defaults to glycine via
#'   [rdkit_mol_from_smiles()].
#' @return Numeric exact mass, or \code{NA_real_}.
#'
#' @examples
#' \dontrun{
#' rdkit_mol_exact_mass()
#' rdkit_mol_exact_mass(rdkit_mol_from_smiles("CCO"))
#' }
#'
#' @export
rdkit_mol_exact_mass <- function(mol = rdkit_mol_from_smiles()) {
  if (is.null(mol)) {
    return(NA_real_)
  }
  Descriptors <- reticulate::import("rdkit.Chem.Descriptors", convert = TRUE)
  as.numeric(Descriptors$ExactMolWt(mol))
}

#' Test whether a mol matches a SMARTS pattern
#'
#' @param mol RDKit molecule. Defaults to glycine via
#'   [rdkit_mol_from_smiles()].
#' @param smarts Character SMARTS pattern. Default \code{"[NH2]"} (primary amine).
#' @param Chem Optional \code{rdkit.Chem} module from [get_RDKit_Chem()].
#' @return Logical scalar.
#'
#' @examples
#' \dontrun{
#' rdkit_has_substruct()
#' rdkit_has_substruct(rdkit_mol_from_smiles("CCO"), "[OH]")
#' }
#'
#' @export
rdkit_has_substruct <- function(mol = rdkit_mol_from_smiles(),
                                smarts = "[NH2]",
                                Chem = get_RDKit_Chem()) {
  if (is.null(mol) || is.null(smarts) || is.na(smarts) || !nzchar(smarts)) {
    return(FALSE)
  }
  pat <- Chem$MolFromSmarts(as.character(smarts))
  if (is.null(pat) || reticulate::py_is_null_xptr(pat)) {
    return(FALSE)
  }
  isTRUE(mol$HasSubstructMatch(pat))
}

#' Plot an RDKit molecule (in-memory Draw)
#'
#' Renders \code{mol} with RDKit \code{Draw.MolToImage} entirely in memory
#' (no temp file) and returns a ggplot that prints like a normal plot.
#'
#' @param mol RDKit molecule. Defaults to glycine via
#'   [rdkit_mol_from_smiles()].
#' @param width,height Integer pixel size passed to RDKit Draw.
#' @return A ggplot object with the molecule raster.
#'
#' @examples
#' \dontrun{
#' rdkit_plot_mol()
#' rdkit_plot_mol(rdkit_mol_from_smiles("CCO"))
#' }
#'
#' @export
rdkit_plot_mol <- function(mol = rdkit_mol_from_smiles(),
                           width = 500,
                           height = 500L) {
  ensure_RDKit_python()
  if (is.null(mol) || reticulate::py_is_null_xptr(mol)) {
    stop("'mol' is NULL or invalid.")
  }
  Draw <- reticulate::import("rdkit.Chem.Draw", convert = TRUE)
  np <- reticulate::import("numpy", convert = TRUE)
  img <- Draw$MolToImage(
    mol,
    size = reticulate::tuple(as.integer(width), as.integer(height))
  )
  arr <- np$asarray(img)
  dims <- dim(arr)
  if (is.null(dims) || length(dims) < 2L) {
    stop("Failed to convert RDKit image to an array.")
  }
  if (length(dims) == 2L) {
    cols <- grDevices::rgb(arr, arr, arr, maxColorValue = 255)
  } else if (dims[3] >= 4L) {
    cols <- grDevices::rgb(arr[, , 1], arr[, , 2], arr[, , 3], arr[, , 4],
                           maxColorValue = 255)
  } else {
    cols <- grDevices::rgb(arr[, , 1], arr[, , 2], arr[, , 3],
                           maxColorValue = 255)
  }
  ras <- grDevices::as.raster(matrix(cols, nrow = dims[1], ncol = dims[2]))
  ggplot2::ggplot() +
    ggplot2::annotation_raster(ras, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void()
}
