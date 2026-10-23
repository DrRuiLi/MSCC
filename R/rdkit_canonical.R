#' RDKit interpreter that matches ICEBERG.
#'
#' @param python Optional path. The default is the ICEBERG CPU environment.
#' @return The python executable.
#' @export
rdkit_canonical_python <- function(python = NULL) {
  if (is.null(python) || !nzchar(python)) {
    python <- Sys.getenv("MSCC_RDKIT_PYTHON", "")
  }
  if (!nzchar(python)) {
    python <- "/mnt/nas/conda/envs/lirui_iceberg_cpu/bin/python"
  }
  if (!file.exists(python)) {
    stop("The ICEBERG RDKit python was not found at ", python, call. = FALSE)
  }
  python
}

.rdkit_canonical_script <- function() {
  path <- system.file("rdkit", "canonical_mol.py", package = "MSCC")
  if (!nzchar(path) || !file.exists(path)) {
    stop("canonical_mol.py is not installed with MSCC.", call. = FALSE)
  }
  path
}

.rdkit_canonical_call <- function(request) {
  script <- .rdkit_canonical_script()
  input <- tempfile("mscc-rdkit-")
  on.exit(unlink(input), add = TRUE)
  jsonlite::write_json(request, input, auto_unbox = TRUE)
  python <- rdkit_canonical_python()
  err <- tempfile("mscc-rdkit-err-")
  on.exit(unlink(err), add = TRUE)
  out <- system2(
    python,
    args = script,
    stdin = input,
    stdout = TRUE,
    stderr = err
  )
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    message <- c(out, readLines(err, warn = FALSE))
    stop(paste(message, collapse = "\n"), call. = FALSE)
  }
  jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE)
}

#' Canonical SMILES from the copied ICEBERG canonical_root.
#'
#' @param smiles One SMILES string.
#' @return The canonical SMILES.
#' @export
rdkit_canonical_smiles <- function(smiles) {
  smiles <- as.character(smiles)[[1]]
  if (!nzchar(smiles) || is.na(smiles)) {
    stop("A SMILES string is required.", call. = FALSE)
  }
  .rdkit_canonical_call(list(task = "canonical", smiles = smiles))$smiles
}

#' Atoms, bonds, formula, and a temporary molblock for one SMILES.
#'
#' The molblock is only for callers that must hand a file to another program.
#' It is not stored on a molecule object.
#'
#' @param smiles One SMILES string.
#' @return A list with smiles, formula, molblock, atoms, and bonds.
#' @export
rdkit_molecule_parts <- function(smiles) {
  smiles <- as.character(smiles)[[1]]
  .rdkit_canonical_call(list(task = "graph", smiles = smiles))
}

#' One maximum common substructure, in canonical atom ids.
#'
#' @param from_smiles Substrate SMILES.
#' @param to_smiles Product SMILES.
#' @return A list with canonical smiles and atom-id pairs.
#' @export
rdkit_mcs_pairs <- function(from_smiles, to_smiles) {
  .rdkit_canonical_call(list(
    task = "mcs",
    from = as.character(from_smiles)[[1]],
    to = as.character(to_smiles)[[1]]
  ))
}

#' RXNMapper atom-map numbers named by canonical atom ids.
#'
#' @param mapped_smiles One mapped SMILES component.
#' @return A list with the canonical smiles and a named list of map numbers.
#' @export
rdkit_rxn_ids <- function(mapped_smiles) {
  .rdkit_canonical_call(list(task = "rxn", smiles = as.character(mapped_smiles)[[1]]))
}

.rdkit_parts_igraph <- function(parts) {
  atoms <- parts$atoms
  if (!length(atoms)) {
    stop("The canonical molecule has no atoms.", call. = FALSE)
  }
  vertices <- data.frame(
    name = vapply(atoms, function(a) a$id, character(1)),
    atom = vapply(atoms, function(a) a$element, character(1)),
    element = vapply(atoms, function(a) a$element, character(1)),
    id = vapply(atoms, function(a) a$id, character(1)),
    label = vapply(atoms, function(a) paste0(" ", a$element, " "), character(1)),
    x = vapply(atoms, function(a) a$x, numeric(1)),
    y = vapply(atoms, function(a) a$y, numeric(1)),
    stringsAsFactors = FALSE
  )
  bonds <- parts$bonds
  if (!length(bonds)) {
    edges <- data.frame(
      from = character(),
      to = character(),
      bond_type = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    edges <- data.frame(
      from = vapply(bonds, function(b) b$from, character(1)),
      to = vapply(bonds, function(b) b$to, character(1)),
      bond_type = vapply(bonds, function(b) b$bond_type, numeric(1)),
      stringsAsFactors = FALSE
    )
  }
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
}

.rdkit_parts_isotopomer <- function(parts) {
  atoms <- parts$atoms
  ids <- vapply(atoms, function(a) a$id, character(1))
  elements <- vapply(atoms, function(a) a$element, character(1))
  row <- as.list(stats::setNames(elements, ids))
  n_carbon <- sum(elements == "C")
  row <- c(
    list(
      isotopomer = paste0(rep("0", n_carbon), collapse = ""),
      isotopologue = "M+0",
      label = "base",
      abundance = 1,
      path = ""
    ),
    row
  )
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
}
