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

#' Map ChemmineR SDF atom IDs to openclatura IUPAC locants
#'
#' Converts an SDF to an in-memory molblock, loads it with RDKit (preserving
#' atomblock order), runs openclatura numbering, and returns a table aligning
#' ChemmineR canonical atom IDs with parent-chain IUPAC locants.
#'
#' Requires a Python env with \code{rdkit} and \code{openclatura} (see
#' [ensure_RDKit_python()]).
#'
#' @section Limitations:
#' Locants come only from openclatura's **parent-skeleton NUMBERING** step
#' (\code{atom_to_locant}). Atoms that are not on that parent (substituents,
#' heteroatoms off the chain, most ring atoms when another fragment is chosen
#' as parent) stay \code{NA}. This is openclatura behavior, not an SDF/RDKit
#' index mismatch.
#'
#' Example: for a steroid acetate named as \emph{\ldots-yl acetate}, openclatura
#' treats the acetate as the parent, so only the two acetate carbons receive
#' locants (\code{C1}/\code{C2}); the tetracyclic carbons remain \code{NA}.
#' Stereo locants embedded in the name string (e.g. \code{1S,2R,11S,\ldots}) are
#' not exported as a per-atom map. Biological numbering (e.g. steroid
#' C1--C17) is a different nomenclature and is not provided here.
#'
#' @param sdf A ChemmineR \code{SDF} or \code{SDFset}.
#' @param condaenv,python Passed to [ensure_RDKit_python()] before import.
#'
#' @return For an \code{SDF}: a data.frame with columns \code{Atom_id},
#'   \code{element}, \code{rdkit_idx}, \code{locant}, \code{IUPAC_id}, and
#'   attribute \code{iupac_name}. For an \code{SDFset}: a named list of such
#'   tables.
#'
#' @examples
#' \dontrun{
#' sdf <- get_smiles_sdf("NCC(O)=O")[[1]]
#' idx <- get_sdf_IUPAC_index(sdf)
#' # C_2 -> C2 (alpha), C_3 -> C1 (carboxyl); N/O often NA
#' attr(idx, "iupac_name")
#' }
#'
#' @export
get_sdf_IUPAC_index <- function(sdf,
                                condaenv = getOption("MSCC.rdkit.condaenv", "env_for_r"),
                                python = getOption("MSCC.rdkit.python", NULL)) {
  if (inherits(sdf, "SDFset")) {
    out <- lapply(seq_along(sdf), function(i) {
      get_sdf_IUPAC_index(sdf[[i]], condaenv = condaenv, python = python)
    })
    nms <- tryCatch(ChemmineR::cid(sdf), error = function(e) NULL)
    if (!is.null(nms) && length(nms) == length(out)) {
      names(out) <- nms
    }
    return(out)
  }
  if (!inherits(sdf, "SDF")) {
    stop("'sdf' must be a ChemmineR SDF or SDFset.")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for get_sdf_IUPAC_index().")
  }
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    stop("Package 'ChemmineR' is required for get_sdf_IUPAC_index().")
  }

  ensure_RDKit_python(condaenv = condaenv, python = python)
  Chem <- get_RDKit_Chem(condaenv = condaenv, python = python)
  oc <- tryCatch(
    reticulate::import("openclatura", convert = TRUE),
    error = function(e) {
      stop(
        "Python package 'openclatura' is required. ",
        "Install with: conda run -n env_for_r python -m pip install openclatura\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  mb <- paste(ChemmineR::sdf2str(sdf), collapse = "\n")
  mol <- Chem$MolFromMolBlock(mb, removeHs = TRUE)
  if (is.null(mol) || reticulate::py_is_null_xptr(mol)) {
    stop("Failed to parse SDF molblock with RDKit MolFromMolBlock().")
  }

  an <- oc$analyze_rdkit_mol(mol)
  decisions <- an$decisions
  num_step <- NULL
  for (step in decisions) {
    phase <- tolower(as.character(step$phase))
    if (length(phase) && any(phase == "numbering")) {
      num_step <- step
      break
    }
  }
  if (is.null(num_step)) {
    stop("openclatura analysis did not include a NUMBERING step.")
  }

  loc_raw <- num_step$data$atom_to_locant
  loc_map <- .openclatura_locant_map(loc_raw)

  ab <- ChemmineR::atomblock(sdf)
  atom_ids <- rownames(ab)
  if (is.null(atom_ids) || !length(atom_ids)) {
    stop("SDF atomblock has no rownames (Atom_id).")
  }
  # ChemmineR bonds()/atomblock element column is typically named via cbind in
  # Molecule helpers; atomblock itself has element in column used by atom labels.
  # Prefer the atom symbol from the bonds table when available.
  bonds_df <- tryCatch(ChemmineR::bonds(sdf), error = function(e) NULL)
  if (!is.null(bonds_df) && "atom" %in% colnames(bonds_df)) {
    elements <- as.character(bonds_df$atom[seq_along(atom_ids)])
  } else {
    # Fallback: parse Element_N style rownames
    elements <- sub("_.*$", "", atom_ids)
  }

  n <- length(atom_ids)
  rdkit_idx <- seq_len(n) - 1L
  locant <- loc_map[as.character(rdkit_idx)]
  locant <- unname(as.character(locant))
  iupac_id <- ifelse(
    is.na(locant) | !nzchar(locant),
    NA_character_,
    paste0(elements, locant)
  )

  out <- data.frame(
    Atom_id = atom_ids,
    element = elements,
    rdkit_idx = rdkit_idx,
    locant = locant,
    IUPAC_id = iupac_id,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "iupac_name") <- as.character(an$name)
  out
}

# Convert openclatura atom_to_locant (Python dict / R list) to a named character vector.
.openclatura_locant_map <- function(loc_raw) {
  if (is.null(loc_raw)) {
    return(setNames(character(), character()))
  }
  if (inherits(loc_raw, "python.builtin.dict") ||
      inherits(loc_raw, "reticulate.python.builtin.dict")) {
    loc_raw <- reticulate::py_to_r(loc_raw)
  }
  if (is.list(loc_raw) || is.vector(loc_raw)) {
    nms <- names(loc_raw)
    if (is.null(nms) || !length(nms)) {
      # integer-keyed list from py_to_r may use [[i]] with names as keys
      nms <- as.character(names(loc_raw))
    }
    vals <- vapply(loc_raw, function(x) {
      if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) {
        return(NA_character_)
      }
      as.character(x)[[1]]
    }, character(1))
    if (is.null(nms) || !any(nzchar(nms))) {
      # keys may be stored only as list names after convert
      nms <- as.character(seq_along(vals) - 1L)
    }
    return(setNames(vals, as.character(nms)))
  }
  setNames(character(), character())
}
