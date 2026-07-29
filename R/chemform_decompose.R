#' MCP-only mass decomposition (neutral exact mass)
#'
#' Enumerate candidate sum formulas that match a **neutral exact mass** within
#' an error window. Uses an imslib/Rdisop-style Money Changing Problem (MCP)
#' solver and intentionally skips isotope distribution calculation/ranking.
#'
#' For ion m/z input with charge, use [chemform_decompose_mz()].
#'
#' @param mass Neutral exact mass (scalar or vector).
#' @param ppm Allowed deviation in ppm.
#' @param mzabs Allowed absolute deviation in Dalton.
#' @param elements Character vector of allowed elements, e.g. `c("C","H","N","O","P","S")`.
#' @param min_elements Minimum element counts. Accepts `NULL` (defaults to 0 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C0H0N0"`.
#' @param max_elements Maximum element counts. Accepts `NULL` (defaults to 999999 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C999H999"`.
#' @param check_rule If `TRUE`, keep only candidates that pass
#'   [chemform_check_seven_golden_rules()] (Rules #1, #2, #4–#6). Default `FALSE`.
#'
#' @return A data.frame with columns `formula`, `exactmass`, `ppm`, and `mass_target`.
#'   Rows are sorted by increasing `abs(ppm)`.
#' @useDynLib MSCC, .registration = TRUE
#' @import Rcpp
#' @export
#' @seealso [chemform_decompose_mz()], [chemform_check_seven_golden_rules()]
chemform_decompose_mass <- function(mass,
                                     ppm = 5,
                                     mzabs = 1e-4,
                                     elements = c("C", "H", "N", "O", "P", "S"),
                                     min_elements = NULL,
                                     max_elements = NULL,
                                     check_rule = T) {

  mass <- as.numeric(mass)
  if (!length(mass)) {
    return(data.frame(
      formula = character(),
      exactmass = numeric(),
      ppm = numeric(),
      mass_target = numeric(),
      stringsAsFactors = FALSE
    ))
  }

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
  max_counts <- as_counts_vec(max_elements, default = 999999L)

  if (any(min_counts > max_counts)) {
    stop("`min_elements` must be <= `max_elements` for each element.")
  }

  call_one <- function(m_target) {
    abs_error <- ppm * m_target * 1e-6 + mzabs
    res <- mcp_decompose_mass(
      mass = m_target,
      abs_error = abs_error,
      mono_masses = as.numeric(mono_masses),
      element_names = elements,
      min_counts = as.integer(min_counts),
      max_counts = as.integer(max_counts)
    )
    if (!length(res$formula)) {
      return(data.frame(
        formula = character(),
        exactmass = numeric(),
        ppm = numeric(),
        mass_target = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    exactmass <- as.numeric(res$exactmass)
    ppm_val <- (exactmass - m_target) / m_target * 1e6
    ord <- order(abs(ppm_val), ppm_val)

    data.frame(
      formula = res$formula[ord],
      exactmass = exactmass[ord],
      ppm = ppm_val[ord],
      mass_target = rep(m_target, length(ord)),
      stringsAsFactors = FALSE
    )
  }

  out_list <- lapply(mass, call_one)
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  if (isTRUE(check_rule) && nrow(out)) {
    keep <- chemform_check_seven_golden_rules(
      chemform = out$formula,
      mass = out$exactmass,
      return = "valid"
    )
    out <- out[keep, , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}


#' MCP mass decomposition from ion m/z
#'
#' Convert ion m/z + charge to neutral exact mass (same electron-mass convention
#' as [chemform_mz()]), then call [chemform_decompose_mass()].
#'
#' Neutral mass conversion when `charge != 0`:
#' `M = mz * abs(charge) + e * charge`, with `e = 0.00054857990943`.
#'
#' @param mz Ion m/z (scalar or vector). If `charge = 0`, treated as neutral mass.
#' @param charge Integer charge `z` (e.g. `+1`, `+2`, `-1`). Use `0` for neutral input.
#' @param ppm Allowed deviation in ppm (applied on the neutral mass used for MCP).
#' @param mzabs Allowed absolute deviation in Dalton.
#' @param elements Character vector of allowed elements.
#' @param min_elements See [chemform_decompose_mass()].
#' @param max_elements See [chemform_decompose_mass()].
#' @param check_rule If `TRUE`, keep only candidates that pass
#'   [chemform_check_seven_golden_rules()] (passed through to
#'   [chemform_decompose_mass()]). Default `FALSE`.
#'
#' @return A data.frame with columns `formula`, `exactmass`, `mz`, `ppm`, `charge`,
#'   and `mz_target`. Rows are sorted by increasing `abs(ppm)` vs the input m/z.
#' @export
#' @seealso [chemform_decompose_mass()], [chemform_mz()], [chemform_check_seven_golden_rules()]
chemform_decompose_mz <- function(mz,
                                   charge = 0,
                                   ppm = 5,
                                   mzabs = 1e-4,
                                   elements = c("C", "H", "N", "O", "P", "S"),
                                   min_elements = NULL,
                                   max_elements = NULL,
                                   check_rule = T) {

  mz <- as.numeric(mz)
  charge <- as.integer(charge)
  if (length(charge) != 1L) stop("`charge` must be a single integer.")

  if (!length(mz)) {
    return(data.frame(
      formula = character(),
      exactmass = numeric(),
      mz = numeric(),
      ppm = numeric(),
      charge = integer(),
      mz_target = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  e_mass <- 0.00054857990943
  abs_charge <- abs(charge)

  mz_to_neutral <- function(m) {
    if (charge == 0L) return(m)
    m * abs_charge + e_mass * charge
  }

  theoretical_mz <- function(exactmass) {
    if (charge == 0L) return(exactmass)
    (exactmass - e_mass * charge) / abs_charge
  }

  call_one <- function(mz_target) {
    M_neutral <- mz_to_neutral(mz_target)
    res <- chemform_decompose_mass(
      mass = M_neutral,
      ppm = ppm,
      mzabs = mzabs,
      elements = elements,
      min_elements = min_elements,
      max_elements = max_elements
    )
    if (!nrow(res)) {
      return(data.frame(
        formula = character(),
        exactmass = numeric(),
        mz = numeric(),
        ppm = numeric(),
        charge = integer(),
        mz_target = numeric(),
        stringsAsFactors = FALSE
      ))
    }

    mz_theory <- theoretical_mz(res$exactmass)
    ppm_val <- (mz_theory - mz_target) / mz_target * 1e6
    ord <- order(abs(ppm_val), ppm_val)

    data.frame(
      formula = res$formula[ord],
      exactmass = res$exactmass[ord],
      mz = mz_theory[ord],
      ppm = ppm_val[ord],
      charge = rep(charge, length(ord)),
      mz_target = rep(mz_target, length(ord)),
      stringsAsFactors = FALSE
    )
  }

  out_list <- lapply(mz, call_one)
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
