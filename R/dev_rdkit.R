#' RDKit bridge via reticulate
#'
#' Helpers for importing RDKit submodules and building molecules from SMILES.
#' Import individual RDKit submodules by name (e.g. \code{get_RDKit_Chem()},
#' or \code{reticulate::import("rdkit.Chem.Descriptors")}) — do not use a
#' single catch-all \code{get_RDKit()} for the whole package.
#'
#' Requires a Python environment with \code{rdkit} configured for
#' \pkg{reticulate}. Call \code{\link{rdkit_reticulate_python}} \emph{before}
#' any other reticulate Python init in the session. Default interpreter:
#' \code{MSCC_RDKIT_RETICULATE_PYTHON}, then \code{MSCC_RDKIT_CLI_PYTHON},
#' then the ICEBERG CPU path (see \code{\link{rdkit_cli_python}}).
#' CLI / subprocess helpers use \code{\link{rdkit_cli_python}} separately.
#'
#' @name dev_rdkit
NULL

#' Bind reticulate to the RDKit Python interpreter (before first import)
#'
#' \code{reticulate::import()} cannot take an env path; the interpreter must
#' be chosen first. On Windows, also adds the conda env
#' \code{Library/bin} directory to the DLL search path so
#' \code{rdMolDraw2D} / Cairo drawing can load. That add is
#' \emph{once per Python process}: each \code{os.add_dll_directory()}
#' cookie is retained, because discarding them leaks Win32 search
#' entries until \code{WinError 206} (the path itself is not too long).
#'
#' Resolution when \code{python} and \code{condaenv} are both unset:
#' \code{MSCC_RDKIT_RETICULATE_PYTHON}, then \code{MSCC_RDKIT_CLI_PYTHON},
#' then the ICEBERG hardcoded path via \code{\link{rdkit_cli_python}}.
#'
#' @param condaenv Optional conda env name or path. Used only when
#'   \code{python} is not set.
#' @param python Optional path to a \code{python} executable. If set, takes
#'   precedence over \code{condaenv} and the env vars above.
#' @return \code{invisible(TRUE)}.
#' @export
rdkit_reticulate_python <- function(condaenv = NULL,
                                    python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for RDKit. Install it and configure a Python env with rdkit.")
  }
  if (!isTRUE(reticulate::py_available(initialize = FALSE))) {
    if (!is.null(python) && nzchar(python)) {
      reticulate::use_python(rdkit_cli_python(python), required = TRUE)
    } else if (!is.null(condaenv) && nzchar(condaenv)) {
      reticulate::use_condaenv(condaenv, required = TRUE)
    } else {
      py <- Sys.getenv("MSCC_RDKIT_RETICULATE_PYTHON", "")
      if (!nzchar(py)) {
        py <- Sys.getenv("MSCC_RDKIT_CLI_PYTHON", "")
      }
      if (nzchar(py)) {
        reticulate::use_python(rdkit_cli_python(py), required = TRUE)
      } else {
        reticulate::use_python(rdkit_cli_python(), required = TRUE)
      }
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
#' @param condaenv,python Passed to [rdkit_reticulate_python()] before import.
#' @return The imported \code{rdkit.Chem} module.
#' @export
get_RDKit_Chem <- function(condaenv = NULL,
                           python = NULL) {
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
  rdkit_reticulate_python(condaenv = condaenv, python = python)
  .rdkit_mod$Chem <- reticulate::py_suppress_warnings(
    reticulate::import("rdkit.Chem", convert = TRUE)
  )
  .rdkit_mod$Chem
}

#' Build an RDKit molecule from SMILES
#'
#' @param smiles Character scalar SMILES. Default \code{"NCC(O)=O"} (glycine),
#'   matching MSCC helpers such as \code{get_Molecule_igraph_from_smiles()}.
#' @param sanitize If \code{TRUE} (default), use RDKit's normal SMILES
#'   parse. If \code{FALSE}, parse without sanitization (needed for some
#'   charged CFM fragment SMILES with non-standard valence) and then run
#'   a soft \code{SanitizeMol(..., catchErrors = TRUE)}.
#' @param Chem Optional \code{rdkit.Chem} module from [get_RDKit_Chem()].
#' @return An RDKit mol object, or \code{NULL} if invalid.
#'
#' @examples
#' \dontrun{
#' mol <- rdkit_mol_from_smiles()
#' mol <- rdkit_mol_from_smiles("CCO")
#' mol <- rdkit_mol_from_smiles("C=CC#CC#[CH2+]", sanitize = FALSE)
#' }
#'
#' @export
rdkit_mol_from_smiles <- function(smiles = "NCC(O)=O",
                                  sanitize = TRUE,
                                  Chem = get_RDKit_Chem()) {
  if (is.null(smiles) || length(smiles) != 1L || is.na(smiles) || !nzchar(smiles)) {
    return(NULL)
  }
  smi <- as.character(smiles)[[1]]
  if (isTRUE(sanitize)) {
    mol <- Chem$MolFromSmiles(smi)
  } else {
    mol <- Chem$MolFromSmiles(smi, sanitize = FALSE)
    if (!is.null(mol) && !reticulate::py_is_null_xptr(mol)) {
      tryCatch(
        Chem$SanitizeMol(mol, catchErrors = TRUE),
        error = function(e) NULL
      )
    }
  }
  # is.null first: py_is_null_xptr(NULL) errors with
  # "Expected a python object, received a NULL".
  if (is.null(mol) || reticulate::py_is_null_xptr(mol)) {
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

#' Plot a Molecule_igraph with RDKit atom ids
#'
#' Draws \code{mol.ig@molecule_info$smiles} with [rdkit_plot_mol()].
#' When \code{show_id} is \code{TRUE}, each atom label is the graph
#' vertex \code{id} (\code{C_1}, \code{N_1}, \ldots). Those ids are the
#' RDKit atom index plus one, prefixed by the element, and they are the
#' column names of \code{@isotopomer}. The label is written as RDKit's
#' \code{atomLabel} property before drawing, which replaces the element
#' symbol.
#'
#' Vertex row order must match the atom order of the canonical SMILES.
#' The function stops if the atom count or element symbols disagree.
#'
#' @param mol.ig A \code{Molecule_igraph}.
#' @param show_id If \code{TRUE}, label atoms with vertex ids. Default
#'   \code{TRUE}.
#' @param highlight_ids Optional atom ids to highlight, or a named
#'   numeric vector of mapping probabilities (names = vertex ids).
#'   Values \code{> 0} are highlighted; colour intensity follows the
#'   probability when names are supplied. Bonds between highlighted
#'   atoms are filled too.
#' @param ... Passed to [rdkit_plot_mol()] / the Cairo drawer. With
#'   \code{show_id = TRUE}, \code{fixedFontSize} defaults to \code{28}
#'   and \code{padding} to \code{0.16} unless those arguments are
#'   supplied.
#' @return A ggplot object with the molecule raster.
#'
#' @examples
#' \dontrun{
#' mig <- get_Molecule_igraph_from_smiles("NCC(O)=O")
#' rdkit_plot_Molecule_igraph(mig)
#' rdkit_plot_Molecule_igraph(mig, highlight_ids = c("C_2", "C_3"))
#' }
#'
#' @export
rdkit_plot_Molecule_igraph <- function(mol.ig,
                                       show_id = TRUE,
                                       highlight_ids = NULL,
                                       ...) {
  if (!inherits(mol.ig, "Molecule_igraph")) {
    stop("'mol.ig' must be a Molecule_igraph.", call. = FALSE)
  }
  smi <- mol.ig@molecule_info$smiles
  if (is.null(smi) || !length(smi) || is.na(smi[[1]]) || !nzchar(smi[[1]])) {
    stop("This molecule has no SMILES to plot.", call. = FALSE)
  }
  mol <- rdkit_mol_from_smiles(smi[[1]])
  if (is.null(mol)) {
    stop("RDKit could not parse the molecule SMILES.", call. = FALSE)
  }
  dots <- list(...)
  if (isTRUE(show_id)) {
    mol <- .rdkit_set_igraph_atom_labels(mol, mol.ig)
    if (is.null(dots$fixedFontSize)) dots$fixedFontSize <- 28
    if (is.null(dots$padding)) dots$padding <- 0.16
  }
  hl <- .rdkit_highlight_from_igraph_ids(mol, mol.ig, highlight_ids)
  draw_args <- list(
    mol = mol,
    width = if (is.null(dots$width)) 800L else dots$width,
    height = if (is.null(dots$height)) 800L else dots$height,
    bondLineWidth = if (is.null(dots$bondLineWidth)) 7 else dots$bondLineWidth,
    fixedFontSize = if (is.null(dots$fixedFontSize)) 80L else dots$fixedFontSize,
    minFontSize = if (is.null(dots$minFontSize)) 72L else dots$minFontSize,
    maxFontSize = if (is.null(dots$maxFontSize)) 88L else dots$maxFontSize,
    padding = if (is.null(dots$padding)) 0.12 else dots$padding,
    additionalAtomLabelPadding = if (is.null(dots$additionalAtomLabelPadding)) {
      0.04
    } else {
      dots$additionalAtomLabelPadding
    },
    multipleBondOffset = if (is.null(dots$multipleBondOffset)) {
      0.16
    } else {
      dots$multipleBondOffset
    },
    scaleBondWidth = isTRUE(dots$scaleBondWidth),
    centreMoleculesBeforeDrawing = if (is.null(dots$centreMoleculesBeforeDrawing)) {
      TRUE
    } else {
      isTRUE(dots$centreMoleculesBeforeDrawing)
    },
    highlight = hl,
    background = if (is.null(dots$background)) c(1, 1, 1) else dots$background
  )
  do.call(.rdkit_draw_mol_ggplot, draw_args)
}

#' Highlight atoms by Molecule_igraph vertex ids / map probabilities.
#'
#' @param mol An RDKit molecule matching \code{mol.ig} atom order.
#' @param mol.ig A \code{Molecule_igraph}.
#' @param highlight_ids Character atom ids, or named numeric probabilities.
#' @return Highlight list for [.rdkit_draw_mol_ggplot()].
#' @keywords internal
.rdkit_highlight_from_igraph_ids <- function(mol, mol.ig, highlight_ids) {
  empty <- list(
    atoms = integer(),
    bonds = integer(),
    atom_colors = reticulate::dict(),
    bond_colors = reticulate::dict()
  )
  if (is.null(highlight_ids) || !length(highlight_ids)) {
    return(empty)
  }
  vd <- vdata(mol.ig)
  if (!"id" %in% colnames(vd)) {
    return(empty)
  }
  graph_ids <- as.character(vd$id)
  if (is.numeric(highlight_ids)) {
    probs <- as.numeric(highlight_ids)
    nms <- names(highlight_ids)
    if (is.null(nms) || !any(nzchar(nms))) {
      return(empty)
    }
    names(probs) <- as.character(nms)
  } else {
    nms <- as.character(highlight_ids)
    nms <- nms[!is.na(nms) & nzchar(nms)]
    if (!length(nms)) return(empty)
    probs <- stats::setNames(rep(1, length(nms)), nms)
  }
  probs[!is.finite(probs) | is.na(probs)] <- 0
  probs <- probs[probs > 0]
  if (!length(probs)) return(empty)

  idx0 <- match(names(probs), graph_ids) - 1L
  keep <- !is.na(idx0)
  idx0 <- as.integer(idx0[keep])
  probs <- probs[keep]
  if (!length(idx0)) return(empty)

  base_rgb <- .rdkit_highlight_palette("co_elute")
  atom_col <- list()
  for (i in seq_along(idx0)) {
    p <- max(0, min(1, as.numeric(probs[[i]])))
    # white -> blue by mapping probability
    rgb <- c(
      1 - (1 - base_rgb[1]) * p,
      1 - (1 - base_rgb[2]) * p,
      1 - (1 - base_rgb[3]) * p
    )
    atom_col[[as.character(idx0[[i]])]] <- rgb
  }

  bond_col <- list()
  if (length(idx0) >= 2L) {
    cmb <- utils::combn(idx0, 2L)
    for (j in seq_len(ncol(cmb))) {
      b <- mol$GetBondBetweenAtoms(
        as.integer(cmb[1L, j]),
        as.integer(cmb[2L, j])
      )
      if (is.null(b) || reticulate::py_is_null_xptr(b)) next
      a1 <- as.character(cmb[1L, j])
      a2 <- as.character(cmb[2L, j])
      rgb <- (atom_col[[a1]] + atom_col[[a2]]) / 2
      bond_col[[as.character(as.integer(b$GetIdx()))]] <- rgb
    }
  }

  atoms <- as.integer(names(atom_col))
  bonds <- as.integer(names(bond_col))
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

.rdkit_set_igraph_atom_labels <- function(mol, mol.ig) {
  vd <- vdata(mol.ig)
  if (!all(c("id", "element") %in% colnames(vd))) {
    stop("Molecule_igraph vertices have no id and element.", call. = FALSE)
  }
  ids <- as.character(vd$id)
  elements <- as.character(vd$element)
  n <- as.integer(mol$GetNumAtoms())
  if (length(ids) != n || length(elements) != n) {
    stop(
      "Molecule_igraph has ", length(ids),
      " atoms and the SMILES has ", n, ".",
      call. = FALSE
    )
  }
  if (any(!nzchar(ids) | is.na(ids))) {
    stop("Molecule_igraph has an empty atom id.", call. = FALSE)
  }
  for (i in seq_len(n)) {
    atom <- mol$GetAtomWithIdx(as.integer(i - 1L))
    symbol <- as.character(atom$GetSymbol())
    if (!identical(symbol, elements[[i]])) {
      stop(
        "Atom ", i, " is ", elements[[i]],
        " on the graph and ", symbol, " in the SMILES.",
        call. = FALSE
      )
    }
    atom$SetProp("atomLabel", ids[[i]])
  }
  mol
}

#' Plot a Molecule_igraph with ggplot geometry
#'
#' Reads 2D coordinates, Kekule bond orders, ring membership, and wedge
#' directions from RDKit after \code{PrepareMolForDrawing()}. Bonds are
#' \code{geom_segment} layers and atom ids are \code{geom_text}, so a PDF
#' device keeps them as vectors. Element colours follow RDKit's default
#' palette. A ring double bond is offset toward the ring centre; a
#' non-ring double bond is a symmetric pair of lines.
#'
#' @param mol.ig A \code{Molecule_igraph}.
#' @param show_id If \code{TRUE}, label every atom with its vertex id.
#'   If \code{FALSE}, label heteroatoms with the element symbol. Default
#'   \code{TRUE}.
#' @param highlight_ids Optional atom ids to highlight, or a named
#'   numeric vector of mapping probabilities (names = vertex ids).
#'   Values \code{> 0} get a filled disc under the atom; colour
#'   intensity follows the probability.
#' @param label_size \code{geom_text} size. Default \code{3.2}.
#' @param bond_linewidth Bond \code{linewidth}. Default \code{0.7}.
#' @param highlight_size Highlight disc size in mm for \code{geom_point}.
#'   Default \code{8}.
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' mig <- get_Molecule_igraph_from_smiles("NCC(O)=O")
#' ggplot_Molecule_igraph(mig)
#' ggplot_Molecule_igraph(mig, highlight_ids = c("C_2", "C_3"))
#' }
#'
#' @export
ggplot_Molecule_igraph <- function(mol.ig,
                                   show_id = TRUE,
                                   highlight_ids = NULL,
                                   label_size = 3.2,
                                   bond_linewidth = 0.7,
                                   highlight_size = 8) {
  drawn <- .ggplot_molecule_drawing(mol.ig, show_id = isTRUE(show_id))
  .ggplot_drawn_plot(
    drawn,
    highlight_ids = highlight_ids,
    label_size = label_size,
    bond_linewidth = bond_linewidth,
    highlight_size = highlight_size
  )
}

#' Plot a SMILES string with ggplot geometry
#'
#' Like [ggplot_Molecule_igraph()], but takes SMILES directly. Useful for
#' CFM fragment SMILES that fail strict sanitization / canonicalization
#' (for example \code{C=CC#CC#[CH2+]}); those are drawn with
#' \code{sanitize = FALSE}.
#'
#' @param smiles Character scalar SMILES.
#' @param atom_label \code{"element"} labels every atom with its element
#'   (and formal charge when non-zero), \code{"id"} uses \code{C_1}-style
#'   ids, \code{"hetero"} hides neutral carbons.
#' @param sanitize Passed to [rdkit_mol_from_smiles()]. If \code{TRUE} and
#'   parsing fails, automatically retries with \code{sanitize = FALSE}.
#' @param highlight_ids Optional atom ids / named probabilities.
#' @param label_size,bond_linewidth,highlight_size Drawing sizes.
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' ggplot_smiles("CCO", atom_label = "element")
#' ggplot_smiles("C=CC#CC#[CH2+]", atom_label = "element")
#' }
#'
#' @export
ggplot_smiles <- function(smiles,
                          atom_label = c("element", "id", "hetero"),
                          sanitize = TRUE,
                          highlight_ids = NULL,
                          label_size = 3.2,
                          bond_linewidth = 0.7,
                          highlight_size = 8) {
  atom_label <- match.arg(atom_label)
  drawn <- .ggplot_smiles_drawing(
    smiles,
    atom_label = atom_label,
    sanitize = isTRUE(sanitize)
  )
  .ggplot_drawn_plot(
    drawn,
    highlight_ids = highlight_ids,
    label_size = label_size,
    bond_linewidth = bond_linewidth,
    highlight_size = highlight_size
  )
}

.ggplot_drawn_plot <- function(drawn,
                               highlight_ids = NULL,
                               label_size = 3.2,
                               bond_linewidth = 0.7,
                               highlight_size = 8) {
  hl <- .ggplot_highlight_atoms(drawn$atoms, highlight_ids, highlight_size)
  p <- ggplot2::ggplot()
  if (nrow(hl)) {
    # Use colour (not fill+colour=NA): shape 21 with NA stroke draws nothing.
    p <- p + ggplot2::geom_point(
      data = hl,
      ggplot2::aes(x = x, y = y, colour = fill, size = size),
      alpha = 0.55
    ) +
      ggplot2::scale_size_identity()
  }
  if (nrow(drawn$segments)) {
    p <- p + ggplot2::geom_segment(
      data = drawn$segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = colour),
      linewidth = bond_linewidth,
      lineend = "butt"
    )
  }
  if (nrow(drawn$wedges)) {
    p <- p + ggplot2::geom_polygon(
      data = drawn$wedges,
      ggplot2::aes(x = x, y = y, group = group, fill = fill),
      colour = NA
    )
  }
  if (nrow(drawn$labels)) {
    p <- p + ggplot2::geom_text(
      data = drawn$labels,
      ggplot2::aes(x = x, y = y, label = label, colour = colour),
      size = label_size
    )
  }
  p +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
}

#' Build highlight points for [ggplot_Molecule_igraph()].
#'
#' @param atoms Atom table from [.ggplot_molecule_drawing()].
#' @param highlight_ids Character ids or named numeric probabilities.
#' @param highlight_size Point size in mm.
#' @return A data.frame with x, y, fill, size (possibly 0 rows).
#' @keywords internal
.ggplot_highlight_atoms <- function(atoms, highlight_ids, highlight_size = 8) {
  empty <- data.frame(
    x = numeric(), y = numeric(), fill = character(), size = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(atoms) || !nrow(atoms) || is.null(highlight_ids) ||
      !length(highlight_ids)) {
    return(empty)
  }
  if (is.numeric(highlight_ids)) {
    probs <- as.numeric(highlight_ids)
    nms <- names(highlight_ids)
    if (is.null(nms) || !any(nzchar(nms))) return(empty)
    names(probs) <- as.character(nms)
  } else {
    nms <- as.character(highlight_ids)
    nms <- nms[!is.na(nms) & nzchar(nms)]
    if (!length(nms)) return(empty)
    probs <- stats::setNames(rep(1, length(nms)), nms)
  }
  probs[!is.finite(probs) | is.na(probs)] <- 0
  probs <- probs[probs > 0]
  if (!length(probs)) return(empty)

  hit <- match(names(probs), as.character(atoms$id))
  keep <- !is.na(hit)
  if (!any(keep)) return(empty)
  hit <- hit[keep]
  probs <- probs[keep]
  base_rgb <- .rdkit_highlight_palette("co_elute")
  fill <- vapply(probs, function(p) {
    p <- max(0, min(1, as.numeric(p)))
    rgb <- c(
      1 - (1 - base_rgb[1]) * p,
      1 - (1 - base_rgb[2]) * p,
      1 - (1 - base_rgb[3]) * p
    )
    grDevices::rgb(rgb[1], rgb[2], rgb[3])
  }, character(1))
  data.frame(
    x = atoms$x[hit],
    y = atoms$y[hit],
    fill = fill,
    size = rep(as.numeric(highlight_size)[1L], length(hit)),
    stringsAsFactors = FALSE
  )
}

.ggplot_element_colour <- function(element) {
  switch(
    as.character(element)[[1]],
    C = "#111111",
    N = "#0000FF",
    O = "#E00000",
    S = "#808000",
    P = "#FF8800",
    F = "#33CCCC",
    Cl = "#00AA00",
    Br = "#992222",
    I = "#660066",
    H = "#888888",
    "#444444"
  )
}

.rdkit_point_xy <- function(conf, idx0) {
  pos <- conf$GetAtomPosition(as.integer(idx0))
  c(as.numeric(pos$x), as.numeric(pos$y))
}

.rdkit_atom_rings <- function(mol) {
  raw <- mol$GetRingInfo()$AtomRings()
  if (is.null(raw) || !length(raw)) {
    return(list())
  }
  if (is.list(raw)) {
    lapply(raw, as.integer)
  } else {
    list(as.integer(raw))
  }
}

.ggplot_atom_labels <- function(ids,
                                elements,
                                charges = NULL,
                                atom_label = c("id", "element", "hetero")) {
  atom_label <- match.arg(atom_label)
  n <- length(elements)
  if (!n) return(character())
  if (identical(atom_label, "id")) {
    return(as.character(ids))
  }
  lab <- as.character(elements)
  if (!is.null(charges) && length(charges) == n) {
    ch <- as.integer(charges)
    ch[!is.finite(ch)] <- 0L
    lab <- ifelse(
      ch == 0L,
      lab,
      ifelse(
        ch == 1L,
        paste0(lab, "+"),
        ifelse(
          ch == -1L,
          paste0(lab, "-"),
          paste0(lab, ifelse(ch > 0L, paste0("+", ch), as.character(ch)))
        )
      )
    )
  }
  if (identical(atom_label, "hetero")) {
    hide <- elements == "C"
    if (!is.null(charges) && length(charges) == n) {
      hide <- hide & (as.integer(charges) == 0L)
    }
    lab[hide] <- ""
  }
  lab
}

.ggplot_mol_drawing <- function(mol,
                                ids,
                                elements,
                                charges = NULL,
                                atom_label = c("id", "element", "hetero")) {
  atom_label <- match.arg(atom_label)
  if (as.integer(mol$GetNumConformers()) < 1L) {
    AllChem <- reticulate::import("rdkit.Chem.AllChem", convert = TRUE)
    AllChem$Compute2DCoords(mol)
  }
  rings <- .rdkit_atom_rings(mol)
  rdMolDraw2D <- reticulate::import(
    "rdkit.Chem.Draw.rdMolDraw2D",
    convert = TRUE
  )
  mol <- tryCatch(
    rdMolDraw2D$PrepareMolForDrawing(mol),
    error = function(e) rdMolDraw2D$PrepareMolForDrawing(mol, kekulize = FALSE)
  )
  conf <- mol$GetConformer()
  n_draw <- as.integer(mol$GetNumAtoms())
  n <- length(elements)
  xy <- vapply(seq_len(n_draw) - 1L, function(i) .rdkit_point_xy(conf, i), numeric(2))
  atoms <- data.frame(
    idx0 = seq_len(n_draw) - 1L,
    id = NA_character_,
    x = xy[1, ],
    y = xy[2, ],
    element = vapply(
      seq_len(n_draw) - 1L,
      function(i) as.character(mol$GetAtomWithIdx(as.integer(i))$GetSymbol()),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  use_n <- min(n, n_draw)
  atoms$id[seq_len(use_n)] <- as.character(ids)[seq_len(use_n)]
  ch <- if (is.null(charges)) {
    NULL
  } else {
    as.integer(charges)[seq_len(use_n)]
  }
  atoms$label <- ""
  atoms$label[seq_len(use_n)] <- .ggplot_atom_labels(
    ids = as.character(ids)[seq_len(use_n)],
    elements = as.character(elements)[seq_len(use_n)],
    charges = ch,
    atom_label = atom_label
  )
  atoms$colour <- vapply(atoms$element, .ggplot_element_colour, character(1))
  bond_len <- .ggplot_bond_lengths(mol, atoms)
  unit_len <- if (length(bond_len)) stats::median(bond_len) else 1.5
  segments <- .ggplot_bond_segments(mol, atoms, rings, unit_len)
  styled <- .ggplot_bond_wedges(mol, atoms, unit_len)
  segments <- rbind(segments, styled$dashes)
  labels <- atoms[nzchar(atoms$label), c("x", "y", "label", "colour"), drop = FALSE]
  list(
    atoms = atoms,
    segments = segments,
    wedges = styled$wedges,
    labels = labels
  )
}

.ggplot_molecule_drawing <- function(mol.ig, show_id) {
  if (!inherits(mol.ig, "Molecule_igraph")) {
    stop("'mol.ig' must be a Molecule_igraph.", call. = FALSE)
  }
  smi <- mol.ig@molecule_info$smiles
  if (is.null(smi) || !length(smi) || is.na(smi[[1]]) || !nzchar(smi[[1]])) {
    stop("This molecule has no SMILES to plot.", call. = FALSE)
  }
  mol <- rdkit_mol_from_smiles(smi[[1]])
  if (is.null(mol)) {
    stop("RDKit could not parse the molecule SMILES.", call. = FALSE)
  }
  vd <- vdata(mol.ig)
  if (!all(c("id", "element") %in% colnames(vd))) {
    stop("Molecule_igraph vertices have no id and element.", call. = FALSE)
  }
  ids <- as.character(vd$id)
  elements <- as.character(vd$element)
  n <- as.integer(mol$GetNumAtoms())
  if (length(ids) != n || length(elements) != n) {
    stop(
      "Molecule_igraph has ", length(ids),
      " atoms and the SMILES has ", n, ".",
      call. = FALSE
    )
  }
  for (i in seq_len(n)) {
    symbol <- as.character(mol$GetAtomWithIdx(as.integer(i - 1L))$GetSymbol())
    if (!identical(symbol, elements[[i]])) {
      stop(
        "Atom ", i, " is ", elements[[i]],
        " on the graph and ", symbol, " in the SMILES.",
        call. = FALSE
      )
    }
  }
  .ggplot_mol_drawing(
    mol = mol,
    ids = ids,
    elements = elements,
    charges = NULL,
    atom_label = if (isTRUE(show_id)) "id" else "hetero"
  )
}

.ggplot_smiles_drawing <- function(smiles,
                                   atom_label = c("element", "id", "hetero"),
                                   sanitize = TRUE) {
  atom_label <- match.arg(atom_label)
  if (is.null(smiles) || length(smiles) != 1L || is.na(smiles) || !nzchar(smiles)) {
    stop("SMILES is empty.", call. = FALSE)
  }
  mol <- rdkit_mol_from_smiles(smiles, sanitize = isTRUE(sanitize))
  if (is.null(mol) && isTRUE(sanitize)) {
    mol <- rdkit_mol_from_smiles(smiles, sanitize = FALSE)
  }
  if (is.null(mol)) {
    stop("RDKit could not parse the SMILES.", call. = FALSE)
  }
  n <- as.integer(mol$GetNumAtoms())
  elements <- character(n)
  charges <- integer(n)
  for (i in seq_len(n)) {
    atom <- mol$GetAtomWithIdx(as.integer(i - 1L))
    elements[[i]] <- as.character(atom$GetSymbol())
    charges[[i]] <- as.integer(atom$GetFormalCharge())
  }
  ids <- paste0(elements, "_", seq_len(n))
  .ggplot_mol_drawing(
    mol = mol,
    ids = ids,
    elements = elements,
    charges = charges,
    atom_label = atom_label
  )
}

.ggplot_bond_lengths <- function(mol, atoms) {
  n_bond <- as.integer(mol$GetNumBonds())
  if (!n_bond) {
    return(numeric())
  }
  vapply(seq_len(n_bond) - 1L, function(i) {
    bond <- mol$GetBondWithIdx(as.integer(i))
    a <- atoms[atoms$idx0 == as.integer(bond$GetBeginAtomIdx()), ]
    b <- atoms[atoms$idx0 == as.integer(bond$GetEndAtomIdx()), ]
    sqrt((b$x - a$x)^2 + (b$y - a$y)^2)
  }, numeric(1))
}

.ggplot_bond_order <- function(bond) {
  ord <- suppressWarnings(as.numeric(bond$GetBondTypeAsDouble()))
  if (!is.finite(ord) || ord < 1.5) {
    1L
  } else if (ord < 2.5) {
    2L
  } else {
    3L
  }
}

.ggplot_bond_style <- function(bond) {
  dir <- toupper(as.character(bond$GetBondDir()))
  if (grepl("WEDGE", dir, fixed = TRUE)) {
    "wedge"
  } else if (grepl("DASH", dir, fixed = TRUE)) {
    "dash"
  } else {
    "plain"
  }
}

.ggplot_shorten <- function(x1, y1, x2, y2, gap1, gap2) {
  dx <- x2 - x1
  dy <- y2 - y1
  len <- sqrt(dx * dx + dy * dy)
  if (!is.finite(len) || len < 1e-8) {
    return(NULL)
  }
  ux <- dx / len
  uy <- dy / len
  if (gap1 + gap2 > len * 0.8) {
    scale <- (len * 0.8) / max(gap1 + gap2, 1e-8)
    gap1 <- gap1 * scale
    gap2 <- gap2 * scale
  }
  list(
    x = x1 + ux * gap1,
    y = y1 + uy * gap1,
    xend = x2 - ux * gap2,
    yend = y2 - uy * gap2,
    px = -uy,
    py = ux,
    len = len
  )
}

.ggplot_ring_centre <- function(rings, i0, j0, atoms) {
  hits <- Filter(function(r) i0 %in% r && j0 %in% r, rings)
  if (!length(hits)) {
    return(NULL)
  }
  ring <- hits[[which.min(lengths(hits))]]
  pts <- atoms[atoms$idx0 %in% ring, c("x", "y"), drop = FALSE]
  if (!nrow(pts)) {
    return(NULL)
  }
  c(mean(pts$x), mean(pts$y))
}

.ggplot_half_segments <- function(x, y, xend, yend, col1, col2) {
  mx <- (x + xend) / 2
  my <- (y + yend) / 2
  data.frame(
    x = c(x, mx),
    y = c(y, my),
    xend = c(mx, xend),
    yend = c(my, yend),
    colour = c(col1, col2),
    stringsAsFactors = FALSE
  )
}

.ggplot_bond_segments <- function(mol, atoms, rings, unit_len) {
  empty <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    colour = character(), stringsAsFactors = FALSE
  )
  n_bond <- as.integer(mol$GetNumBonds())
  if (!n_bond) {
    return(empty)
  }
  gap_for <- function(idx0) {
    lab <- atoms$label[atoms$idx0 == idx0]
    if (length(lab) && nzchar(lab[[1]])) 0.22 * unit_len else 0
  }
  pieces <- list()
  for (i in seq_len(n_bond) - 1L) {
    bond <- mol$GetBondWithIdx(as.integer(i))
    if (.ggplot_bond_style(bond) != "plain") {
      next
    }
    i0 <- as.integer(bond$GetBeginAtomIdx())
    j0 <- as.integer(bond$GetEndAtomIdx())
    a <- atoms[atoms$idx0 == i0, ]
    b <- atoms[atoms$idx0 == j0, ]
    order <- .ggplot_bond_order(bond)
    base <- .ggplot_shorten(a$x, a$y, b$x, b$y, gap_for(i0), gap_for(j0))
    if (is.null(base)) {
      next
    }
    centre <- .ggplot_ring_centre(rings, i0, j0, atoms)
    lines <- .ggplot_offset_lines(base, order, centre, unit_len)
    for (k in seq_len(nrow(lines))) {
      pieces[[length(pieces) + 1L]] <- .ggplot_half_segments(
        lines$x[k], lines$y[k], lines$xend[k], lines$yend[k],
        a$colour, b$colour
      )
    }
  }
  if (!length(pieces)) {
    return(empty)
  }
  do.call(rbind, pieces)
}

.ggplot_offset_lines <- function(base, order, centre, unit_len) {
  one <- data.frame(
    x = base$x, y = base$y, xend = base$xend, yend = base$yend,
    stringsAsFactors = FALSE
  )
  if (order < 2L) {
    return(one)
  }
  sep <- 0.16 * unit_len
  shift <- function(line, dist) {
    data.frame(
      x = line$x + base$px * dist,
      y = line$y + base$py * dist,
      xend = line$xend + base$px * dist,
      yend = line$yend + base$py * dist,
      stringsAsFactors = FALSE
    )
  }
  if (order >= 3L) {
    return(rbind(shift(one, sep), one, shift(one, -sep)))
  }
  if (is.null(centre)) {
    return(rbind(shift(one, sep / 2), shift(one, -sep / 2)))
  }
  mid_x <- (base$x + base$xend) / 2
  mid_y <- (base$y + base$yend) / 2
  sign <- if ((centre[1] - mid_x) * base$px + (centre[2] - mid_y) * base$py >= 0) {
    1
  } else {
    -1
  }
  inner_gap <- 0.14 * unit_len
  dx <- base$xend - base$x
  dy <- base$yend - base$y
  len <- sqrt(dx * dx + dy * dy)
  if (len > 1e-8 && 2 * inner_gap < len * 0.7) {
    ux <- dx / len
    uy <- dy / len
    inner <- data.frame(
      x = base$x + ux * inner_gap,
      y = base$y + uy * inner_gap,
      xend = base$xend - ux * inner_gap,
      yend = base$yend - uy * inner_gap,
      stringsAsFactors = FALSE
    )
  } else {
    inner <- one
  }
  rbind(one, shift(inner, sign * sep))
}

.ggplot_bond_wedges <- function(mol, atoms, unit_len) {
  empty_w <- data.frame(
    x = numeric(), y = numeric(), group = integer(), fill = character(),
    stringsAsFactors = FALSE
  )
  empty_d <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    colour = character(), stringsAsFactors = FALSE
  )
  n_bond <- as.integer(mol$GetNumBonds())
  wedges <- list()
  dashes <- list()
  for (i in seq_len(n_bond) - 1L) {
    bond <- mol$GetBondWithIdx(as.integer(i))
    style <- .ggplot_bond_style(bond)
    if (style == "plain") {
      next
    }
    i0 <- as.integer(bond$GetBeginAtomIdx())
    j0 <- as.integer(bond$GetEndAtomIdx())
    a <- atoms[atoms$idx0 == i0, ]
    b <- atoms[atoms$idx0 == j0, ]
    gap1 <- if (nzchar(a$label)) 0.22 * unit_len else 0
    gap2 <- if (nzchar(b$label)) 0.22 * unit_len else 0.08 * unit_len
    base <- .ggplot_shorten(a$x, a$y, b$x, b$y, gap1, gap2)
    if (is.null(base)) {
      next
    }
    if (style == "wedge") {
      wide <- 0.18 * unit_len
      wedges[[length(wedges) + 1L]] <- data.frame(
        x = c(base$x, base$xend + base$px * wide, base$xend - base$px * wide),
        y = c(base$y, base$yend + base$py * wide, base$yend - base$py * wide),
        group = length(wedges) + 1L,
        fill = a$colour,
        stringsAsFactors = FALSE
      )
    } else {
      n_dash <- 6L
      dx <- base$xend - base$x
      dy <- base$yend - base$y
      for (k in seq_len(n_dash)) {
        if (k %% 2L == 0L) {
          next
        }
        t0 <- (k - 1) / n_dash
        t1 <- (k - 0.25) / n_dash
        dashes[[length(dashes) + 1L]] <- .ggplot_half_segments(
          base$x + dx * t0, base$y + dy * t0,
          base$x + dx * t1, base$y + dy * t1,
          a$colour, b$colour
        )
      }
    }
  }
  list(
    wedges = if (length(wedges)) do.call(rbind, wedges) else empty_w,
    dashes = if (length(dashes)) do.call(rbind, dashes) else empty_d
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
  rdkit_reticulate_python()
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
  rdkit_reticulate_python()
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
#' [rdkit_reticulate_python()]).
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
#' @param condaenv,python Passed to [rdkit_reticulate_python()] before import.
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
                                condaenv = NULL,
                                python = NULL) {
  if (length(sdf) > 1L && is.character(sdf)) {
    return(lapply(sdf, get_sdf_IUPAC_index, condaenv = condaenv, python = python))
  }
  if (!is.character(sdf) || length(sdf) != 1L || !nzchar(sdf)) {
    stop("Pass a SMILES string to get_sdf_IUPAC_index().", call. = FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for get_sdf_IUPAC_index().")
  }

  parts <- rdkit_molecule_parts(sdf)
  rdkit_reticulate_python(condaenv = condaenv, python = python)
  Chem <- get_RDKit_Chem(condaenv = condaenv, python = python)
  oc <- tryCatch(
    reticulate::import("openclatura", convert = TRUE),
    error = function(e) {
      stop(
        "Python package 'openclatura' is required. ",
        "Install with: ", rdkit_cli_python(),
        " -m pip install openclatura\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  mol <- Chem$MolFromSmiles(parts$smiles)
  if (is.null(mol) || reticulate::py_is_null_xptr(mol)) {
    stop("Failed to parse the canonical SMILES with RDKit MolFromSmiles().")
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

  atom_ids <- vapply(parts$atoms, function(a) a$id, character(1))
  elements <- vapply(parts$atoms, function(a) a$element, character(1))
  if (!length(atom_ids)) {
    stop("The canonical molecule has no atoms.", call. = FALSE)
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
