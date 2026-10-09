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
#' \code{rdMolDraw2D} / Cairo drawing can load. That add is
#' \emph{once per Python process}: each \code{os.add_dll_directory()}
#' cookie is retained, because discarding them leaks Win32 search
#' entries until \code{WinError 206} (the path itself is not too long).
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
# Keep the AddDllDirectory cookies on sys: the returned objects close on
# GC, and repeating add_dll_directory (~4095 times) raises WinError 206
# even for a short path such as .../Library/bin.
.rdkit_add_windows_dll_dirs <- function() {
  if (.Platform$OS.type != "windows") {
    return(invisible(FALSE))
  }
  if (!isTRUE(reticulate::py_available(initialize = TRUE))) {
    return(invisible(FALSE))
  }
  tryCatch(
    reticulate::py_run_string(
      paste(
        "import os",
        "from pathlib import Path",
        "import sys",
        "if not hasattr(sys, '_mscc_rdkit_dll_cookies'):",
        "    sys._mscc_rdkit_dll_cookies = {}",
        "for p in (Path(sys.prefix) / 'Library' / 'bin', Path(sys.base_prefix) / 'Library' / 'bin'):",
        "    if not p.is_dir():",
        "        continue",
        "    key = os.path.normcase(os.path.abspath(str(p)))",
        "    if key in sys._mscc_rdkit_dll_cookies:",
        "        continue",
        "    try:",
        "        sys._mscc_rdkit_dll_cookies[key] = os.add_dll_directory(str(p))",
        "    except OSError:",
        "        prepend = str(p)",
        "        path = os.environ.get('PATH', '')",
        "        if prepend.lower() not in path.lower():",
        "            os.environ['PATH'] = prepend + os.pathsep + path",
        sep = "\n"
      )
    ),
    error = function(e) invisible(NULL)
  )
  invisible(TRUE)
}

.rdkit_mod <- new.env(parent = emptyenv())

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
  chem <- .rdkit_mod$Chem
  if (!is.null(chem)) {
    dead <- tryCatch(
      reticulate::py_is_null_xptr(chem),
      error = function(e) TRUE
    )
    if (!isTRUE(dead)) {
      return(chem)
    }
  }
  ensure_RDKit_python(condaenv = condaenv, python = python)
  .rdkit_mod$Chem <- reticulate::py_suppress_warnings(
    reticulate::import("rdkit.Chem", convert = TRUE)
  )
  .rdkit_mod$Chem
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

#' Plot an RDKit molecule (in-memory Cairo Draw)
#'
#' Renders \code{mol} with RDKit \code{MolDraw2DCairo} entirely in memory
#' (no temp file) and returns a ggplot that prints like a normal plot.
#' Defaults use a square canvas, thick bonds, and large C/O/N labels.
#'
#' @param mol RDKit molecule. Defaults to glycine via
#'   [rdkit_mol_from_smiles()].
#' @param width,height Integer pixel size of the Cairo canvas. Default
#'   \code{800}.
#' @param bondLineWidth Bond pen width. Default \code{7}.
#' @param fixedFontSize Atom-label size in pixels. Use \code{-1} to let
#'   RDKit scale with the molecule. Default \code{80}.
#' @param minFontSize,maxFontSize Font-size bounds when
#'   \code{fixedFontSize} is unused. Defaults \code{72} and \code{88}.
#' @param padding Fraction of the canvas left as margin. Default
#'   \code{0.12}.
#' @param additionalAtomLabelPadding Extra gap around heteroatom labels.
#'   Default \code{0.04}.
#' @param multipleBondOffset Offset between lines of a double/triple
#'   bond. Default \code{0.16}.
#' @param scaleBondWidth If \code{TRUE}, bond width scales with molecule
#'   size. Default \code{FALSE}.
#' @param centreMoleculesBeforeDrawing Centre the drawing in the canvas.
#'   Default \code{TRUE}.
#' @return A ggplot object with the molecule raster.
#'
#' @examples
#' \dontrun{
#' rdkit_plot_mol()
#' rdkit_plot_mol(rdkit_mol_from_smiles("CCO"))
#' rdkit_plot_mol(rdkit_mol_from_smiles("CCO"), fixedFontSize = 48)
#' }
#'
#' @export
rdkit_plot_mol <- function(mol = rdkit_mol_from_smiles(),
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
                           centreMoleculesBeforeDrawing = TRUE) {
  .rdkit_draw_mol_ggplot(
    mol = mol,
    width = width,
    height = height,
    bondLineWidth = bondLineWidth,
    fixedFontSize = fixedFontSize,
    minFontSize = minFontSize,
    maxFontSize = maxFontSize,
    padding = padding,
    additionalAtomLabelPadding = additionalAtomLabelPadding,
    multipleBondOffset = multipleBondOffset,
    scaleBondWidth = scaleBondWidth,
    centreMoleculesBeforeDrawing = centreMoleculesBeforeDrawing,
    highlight = NULL
  )
}

#' Plot an RDKit molecule with SMARTS highlights
#'
#' Same Cairo drawing as [rdkit_plot_mol()], but fills atoms and bonds
#' that match \code{highlight_smarts}. Named roles
#' \code{structure} / \code{co_elute} / \code{product} use gold, blue,
#' and green.
#'
#' @inheritParams rdkit_plot_mol
#' @param highlight_smarts SMARTS string or character vector. Matching
#'   atoms and bonds are highlighted. Named values
#'   \code{structure} / \code{co_elute} / \code{product} select the
#'   palette. Unnamed values use gold.
#' @param background RGB vector in \code{[0, 1]}. Length 3 or 4
#'   (alpha). Default white.
#' @return A ggplot object with the molecule raster.
#'
#' @examples
#' \dontrun{
#' mol <- rdkit_mol_from_smiles("C(=O)(N)N")
#' rdkit_plot_mol_highlight(mol, c(structure = "[NX3H2][CX3]=O"))
#' }
#'
#' @export
rdkit_plot_mol_highlight <- function(mol = rdkit_mol_from_smiles(),
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
                                     background = c(1, 1, 1)) {
  ensure_RDKit_python()
  hl <- .rdkit_highlight_from_smarts(mol, highlight_smarts)
  .rdkit_draw_mol_ggplot(
    mol = mol,
    width = width,
    height = height,
    bondLineWidth = bondLineWidth,
    fixedFontSize = fixedFontSize,
    minFontSize = minFontSize,
    maxFontSize = maxFontSize,
    padding = padding,
    additionalAtomLabelPadding = additionalAtomLabelPadding,
    multipleBondOffset = multipleBondOffset,
    scaleBondWidth = scaleBondWidth,
    centreMoleculesBeforeDrawing = centreMoleculesBeforeDrawing,
    highlight = hl,
    background = background
  )
}

.rdkit_draw_mol_ggplot <- function(mol,
                                   width,
                                   height,
                                   bondLineWidth,
                                   fixedFontSize,
                                   minFontSize,
                                   maxFontSize,
                                   padding,
                                   additionalAtomLabelPadding,
                                   multipleBondOffset,
                                   scaleBondWidth,
                                   centreMoleculesBeforeDrawing,
                                   highlight = NULL,
                                   background = c(1, 1, 1)) {
  ensure_RDKit_python()
  if (is.null(mol) || reticulate::py_is_null_xptr(mol)) {
    stop("'mol' is NULL or invalid.")
  }
  width <- as.integer(width)[1L]
  height <- as.integer(height)[1L]
  if (!is.finite(width) || !is.finite(height) || width < 1L || height < 1L) {
    stop("'width' and 'height' must be positive integers.")
  }
  AllChem <- reticulate::import("rdkit.Chem.AllChem", convert = TRUE)
  rdMolDraw2D <- reticulate::import(
    "rdkit.Chem.Draw.rdMolDraw2D",
    convert = TRUE
  )
  np <- reticulate::import("numpy", convert = TRUE)
  io <- reticulate::import("io", convert = TRUE)
  Image <- reticulate::import("PIL.Image", convert = TRUE)
  if (as.integer(mol$GetNumConformers()) < 1L) {
    AllChem$Compute2DCoords(mol)
  }
  mol <- rdMolDraw2D$PrepareMolForDrawing(mol)
  drawer <- rdMolDraw2D$MolDraw2DCairo(width, height)
  opts <- drawer$drawOptions()
  .rdkit_set_draw_opt(opts, "bondLineWidth", as.numeric(bondLineWidth)[1L])
  .rdkit_set_draw_opt(opts, "fixedFontSize", as.integer(fixedFontSize)[1L])
  .rdkit_set_draw_opt(opts, "minFontSize", as.integer(minFontSize)[1L])
  .rdkit_set_draw_opt(opts, "maxFontSize", as.integer(maxFontSize)[1L])
  .rdkit_set_draw_opt(opts, "padding", as.numeric(padding)[1L])
  .rdkit_set_draw_opt(
    opts, "additionalAtomLabelPadding",
    as.numeric(additionalAtomLabelPadding)[1L]
  )
  .rdkit_set_draw_opt(
    opts, "multipleBondOffset",
    as.numeric(multipleBondOffset)[1L]
  )
  .rdkit_set_draw_opt(
    opts, "scaleBondWidth",
    isTRUE(scaleBondWidth)
  )
  .rdkit_set_draw_opt(
    opts, "centreMoleculesBeforeDrawing",
    isTRUE(centreMoleculesBeforeDrawing)
  )
  has_hl <- is.list(highlight) && length(highlight$atoms)
  .rdkit_set_draw_opt(opts, "fillHighlights", isTRUE(has_hl))
  if (has_hl) {
    .rdkit_set_draw_opt(opts, "highlightRadius", 0.45)
    .rdkit_set_draw_opt(opts, "highlightBondWidthMultiplier", 16)
  }
  .rdkit_set_background(opts, background)
  if (has_hl) {
    drawer$DrawMolecule(
      mol,
      highlightAtoms = as.integer(highlight$atoms),
      highlightBonds = as.integer(highlight$bonds),
      highlightAtomColors = highlight$atom_colors,
      highlightBondColors = highlight$bond_colors
    )
  } else {
    drawer$DrawMolecule(mol)
  }
  drawer$FinishDrawing()
  img <- Image$open(io$BytesIO(drawer$GetDrawingText()))
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
    ggplot2::coord_fixed(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      ratio = height / width
    ) +
    ggplot2::theme_void()
}

.rdkit_highlight_palette <- function(name) {
  key <- tolower(gsub("_SMARTS$", "", as.character(name)[1L], perl = TRUE))
  switch(
    key,
    "structure" = c(1.00, 0.878, 0.541),
    "co_elute" = c(0.565, 0.792, 0.976),
    "product" = c(0.647, 0.839, 0.655),
    c(1.00, 0.878, 0.541)
  )
}

.rdkit_highlight_from_smarts <- function(mol, smarts) {
  empty <- list(
    atoms = integer(),
    bonds = integer(),
    atom_colors = reticulate::dict(),
    bond_colors = reticulate::dict()
  )
  if (is.null(smarts)) {
    return(empty)
  }
  nms <- names(smarts)
  smarts <- as.character(smarts)
  if (is.null(nms)) {
    nms <- rep("structure", length(smarts))
  }
  keep <- !is.na(smarts) & nzchar(trimws(smarts))
  smarts <- smarts[keep]
  nms <- nms[keep]
  if (!length(smarts)) {
    return(empty)
  }
  Chem <- get_RDKit_Chem()
  atom_col <- list()
  bond_col <- list()
  for (i in seq_along(smarts)) {
    pat <- Chem$MolFromSmarts(smarts[[i]])
    if (is.null(pat) || reticulate::py_is_null_xptr(pat)) {
      next
    }
    hits <- mol$GetSubstructMatches(pat)
    hits <- tryCatch(reticulate::py_to_r(hits), error = function(e) hits)
    if (!length(hits)) {
      next
    }
    rgb <- .rdkit_highlight_palette(nms[[i]])
    if (!is.list(hits)) {
      hits <- list(hits)
    }
    for (h in hits) {
      h <- as.integer(unlist(h, use.names = FALSE))
      for (a in h) {
        atom_col[[as.character(a)]] <- rgb
      }
      if (length(h) >= 2L) {
        cmb <- utils::combn(h, 2L)
        for (j in seq_len(ncol(cmb))) {
          b <- mol$GetBondBetweenAtoms(
            as.integer(cmb[1L, j]),
            as.integer(cmb[2L, j])
          )
          if (!is.null(b) && !reticulate::py_is_null_xptr(b)) {
            bond_col[[as.character(as.integer(b$GetIdx()))]] <- rgb
          }
        }
      }
    }
  }
  atoms <- as.integer(names(atom_col))
  bonds <- as.integer(names(bond_col))
  if (!length(atoms)) {
    return(empty)
  }
  atom_colors <- reticulate::dict()
  for (nm in names(atom_col)) {
    rgb <- atom_col[[nm]]
    atom_colors[as.integer(nm)] <- reticulate::tuple(rgb[1], rgb[2], rgb[3])
  }
  bond_colors <- reticulate::dict()
  for (nm in names(bond_col)) {
    rgb <- bond_col[[nm]]
    bond_colors[as.integer(nm)] <- reticulate::tuple(rgb[1], rgb[2], rgb[3])
  }
  list(
    atoms = atoms,
    bonds = bonds,
    atom_colors = atom_colors,
    bond_colors = bond_colors
  )
}

.rdkit_set_background <- function(opts, rgb) {
  rgb <- as.numeric(rgb)
  if (!length(rgb) || any(!is.finite(rgb))) {
    rgb <- c(1, 1, 1)
  }
  r <- rgb[1L]
  g <- if (length(rgb) >= 2L) rgb[2L] else r
  b <- if (length(rgb) >= 3L) rgb[3L] else r
  a <- if (length(rgb) >= 4L) rgb[4L] else 1
  col <- reticulate::tuple(
    as.numeric(r), as.numeric(g), as.numeric(b), as.numeric(a)
  )
  ok <- tryCatch({
    opts$setBackgroundColour(col)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    ok <- .rdkit_set_draw_opt(opts, "bgColour", col)
  }
  ok
}

.rdkit_set_draw_opt <- function(opts, name, value) {
  tryCatch(
    {
      opts[[name]] <- value
      TRUE
    },
    error = function(e) FALSE
  )
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
