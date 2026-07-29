#' MCP-only mass decomposition (formula calculator)
#'
#' Enumerate candidate sum formulas that match a neutral exact mass (or an ion m/z
#' converted to neutral mass) within an error window. This implementation uses
#' an imslib/Rdisop-style Money Changing Problem (MCP) solver, and intentionally
#' skips isotope distribution calculation/ranking.
#'
#' @param mass Exact neutral mass (`charge = 0`) or ion m/z (`charge != 0`).
#' @param ppm Allowed deviation in ppm.
#' @param mzabs Allowed absolute deviation in Dalton.
#' @param charge Integer charge `z` (e.g. `+1`, `+2`, `-1`). Use `0` for neutral mass input.
#' @param elements Character vector of allowed elements, e.g. `c("C","H","N","O","P","S")`.
#' @param min_elements Minimum element counts. Accepts `NULL` (defaults to 0 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C0H0N0"`.
#' @param max_elements Maximum element counts. Accepts `NULL` (defaults to 999 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C999H999"`.
#'
#' @return A data.frame with columns `formula`, `exactmass`, `mz`, `ppm`, and `charge`.
#'   Rows are sorted by increasing `abs(ppm)`.
#' @useDynLib MSCC, .registration = TRUE
#' @import Rcpp
#' @export
chemform_decompose_mass <- function(mass,
                                     ppm = 5,
                                     mzabs = 1e-4,
                                     charge = 0,
                                     elements = c("C", "H", "N", "O", "P", "S"),
                                     min_elements = NULL,
                                     max_elements = NULL) {

  mass <- as.numeric(mass)
  if (!length(mass)) {
    return(data.frame(formula = character(), exactmass = numeric(),
                      mz = numeric(), ppm = numeric(), charge = integer(),
                      mass_target = numeric(), stringsAsFactors = FALSE))
  }
  charge <- as.integer(charge)
  if (length(charge) != 1L) stop("`charge` must be a single integer.")

  # Canonicalize elements (drop isotope prefixes if any).
  elements <- vapply(elements, MSCC::get_ele_uniso, character(1))

  # Harmonize elem table schema and pull mono masses.
  et <- MSCC::get_elem_table()
  if (!all(c("symbol", "mass", "is.isotope") %in% colnames(et))) {
    stop("`MSCC::elem_table` is missing expected columns: symbol, mass, is.isotope.")
  }

  et_uni <- et[!as.logical(et$is.isotope), , drop = FALSE]
  mono_masses <- et_uni$mass[match(elements, as.character(et_uni$symbol))]

  if (any(is.na(mono_masses))) {
    stop("Some `elements` were not found in monoisotopic element table.")
  }

  n <- length(elements)

  # Translate min/max constraints to integer vectors aligned with `elements`.
  as_counts_vec <- function(x, default) {
    if (is.null(x)) return(rep(default, n))
    if (is.character(x) && length(x) == 1L) {
      # parse formula and extract counts per element
      mat <- chemform_parse(x, return = "matrix")
      counts <- vapply(elements, function(e) {
        if (e %in% colnames(mat)) as.integer(mat[, e, drop = TRUE]) else 0L
      }, integer(1))
      return(as.integer(counts))
    }
    if (is.numeric(x) && !is.null(names(x))) {
      v <- as.integer(x)
      out <- rep(default, n)
      idx <- match(elements, names(v))
      out[!is.na(idx)] <- v[idx[!is.na(idx)]]
      return(out)
    }
    if (is.numeric(x) && length(x) == n) {
      return(as.integer(x))
    }
    stop("`min_elements`/`max_elements` must be NULL, named integer vector, or a single formula string.")
  }

  min_counts <- as_counts_vec(min_elements, default = 0L)
  max_default <- 999999L
  max_counts <- as_counts_vec(max_elements, default = max_default)

  if (any(min_counts > max_counts)) {
    stop("`min_elements` must be <= `max_elements` for each element.")
  }

  # MCP error window is expressed on neutral exact mass.
  e_mass <- 0.00054857990943
  abs_charge <- abs(charge)

  convert_to_neutral <- function(m) {
    if (charge == 0L) return(m)
    m * abs_charge + e_mass * charge
  }

  to_theoretical_mz <- function(exactmass) {
    if (charge == 0L) return(exactmass)
    (exactmass - e_mass * charge) / abs_charge
  }

  # Call the compiled MCP-only solver.
  call_one <- function(m_target) {
    M_neutral <- convert_to_neutral(m_target)
    abs_error <- ppm * M_neutral * 1e-6 + mzabs
    res <- mcp_decompose_mass(
      mass = M_neutral,
      abs_error = abs_error,
      mono_masses = as.numeric(mono_masses),
      element_names = elements,
      min_counts = as.integer(min_counts),
      max_counts = as.integer(max_counts)
    )
    if (!length(res$formula)) {
      return(data.frame(formula = character(),
                         exactmass = numeric(),
                         mz = numeric(),
                         ppm = numeric(),
                         charge = integer(),
                         mass_target = numeric(),
                         stringsAsFactors = FALSE))
    }
    exactmass <- as.numeric(res$exactmass)
    mz_theory <- to_theoretical_mz(exactmass)

    ppm_val <- (mz_theory - m_target) / m_target * 1e6
    ord <- order(abs(ppm_val), ppm_val)

    data.frame(
      formula = res$formula[ord],
      exactmass = exactmass[ord],
      mz = mz_theory[ord],
      ppm = ppm_val[ord],
      charge = rep(charge, length(ord)),
      mass_target = rep(m_target, length(ord)),
      stringsAsFactors = FALSE
    )
  }

  out_list <- lapply(mass, call_one)
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

