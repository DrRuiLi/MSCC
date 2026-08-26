#' MCP-only mass decomposition (neutral exact mass)
#'
#' Enumerate candidate sum formulas that match a **neutral exact mass** within
#' an error window. Uses an imslib/Rdisop-style Money Changing Problem (MCP)
#' solver and intentionally skips isotope distribution calculation/ranking.
#'
#' For ion m/z input with charge, use [chemform_decompose_mz()].
#'
#' @param mass Neutral exact mass (scalar or vector).
#' @param ppm Allowed deviation in ppm. `NULL` is the same as `Inf` (unused).
#'   Default `5`. Combined with `mzabs` by taking the **tighter** (smaller)
#'   window; they are not added. See **Mass tolerance**.
#' @param mzabs Allowed absolute deviation in Dalton. `NULL` is the same as
#'   `Inf` (unused). Default `NULL` (ppm-only). Combined with `ppm` by taking
#'   the tighter window (see **Mass tolerance**).
#' @param elements Character vector of allowed elements, e.g. `c("C","H","N","O","P","S")`.
#' @param min_elements Minimum element counts. Accepts `NULL` (defaults to 0 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C0H0N0"`.
#'   Counts may be **negative** (e.g. `c(H = -3)` or `"H-3"`) for replacement / difference
#'   formulas such as COOH → COONa (`H-1Na`). Internally the target mass is shifted so the
#'   MCP solver still enumerates non-negative relative counts.
#' @param max_elements Maximum element counts. Accepts `NULL` (defaults to 999999 for each element),
#'   a named integer vector (names are element symbols), or a single formula string like `"C999H999"`.
#'   Applied as a **post-filter** after MCP enumeration (does not prune the search). With the
#'   default, mass matching already limits atom counts, so this is usually a no-op unless you
#'   set tighter chemical priors (e.g. `c(C = 10, O = 5)`).
#' @param check_rule If `TRUE`, keep only candidates that pass
#'   [chemform_check_seven_golden_rules()] (Rules #1, #2, #4–#6). Default `FALSE`.
#'   Automatically disabled (with a message) when any `min_elements` count is negative,
#'   because the golden rules assume molecular formulas, not signed replacements.
#'
#' @details
#' ## Mass tolerance
#' The C++ MCP solver takes a single absolute window `abs_error` (Dalton).
#' `NULL` and `Inf` are equivalent: that constraint is unused. `ppm` is
#' converted to Dalton in R, then the **minimum** of the two windows is
#' passed through (unused arguments treated as `Inf`):
#'
#' `abs_error = min(ppm * |mass| * 1e-6, mzabs)`
#'
#' Defaults are `ppm = 5` and `mzabs = NULL`, so the default search is 5 ppm.
#' To use only an absolute window, set `ppm = NULL` (or `Inf`). If both are
#' finite, the tighter window is used. A finite `0` collapses that side to
#' exact match (within solver precision). At least one of `ppm` or `mzabs`
#' must be finite.
#'
#' ## Enumeration and cost
#' With non-negative `min_elements` (the default), MCP enumerates **all** compositions of the
#' allowed `elements` whose exact mass falls in the tolerance window. There is no separate
#' `floor(mass / mono_mass)` cap in the C++ core: remaining mass during recursion already
#' bounds each count (e.g. mass 13 cannot yield two carbons). Cost tracks the **number of
#' solutions** (and alphabet size), not the default `max_elements = 999999`. Typical CHNO /
#' CHNOPS masses are inexpensive; large mass, wide `ppm`/`mzabs`, many light atoms (especially H),
#' or a large alphabet increase cost because more formulas fit the window.
#'
#' @return A data.frame with columns `formula`, `exactmass`, `ppm`, and `mass_target`.
#'   Rows are sorted by increasing `abs(ppm)`.
#'
#' @examples
#' # Neutral molecule (non-negative counts)
#' chemform_decompose_mass(180.0634, ppm = 5, check_rule = FALSE)
#'
#' # Replacement delta: COOH -> COONa is H-1Na (mass change Na - H)
#' dM <- chemform_mz("Na") - chemform_mz("H")
#' chemform_decompose_mass(
#'   mass = dM,
#'   elements = c("H", "C", "O", "Na"),
#'   min_elements = c(H = -3),
#'   max_elements = c(H = 10, Na = 3),
#'   check_rule = FALSE
#' )
#'
#' @useDynLib MSCC, .registration = TRUE
#' @import Rcpp
#' @export
#' @seealso [chemform_decompose_mz()], [chemform_check_seven_golden_rules()]
chemform_decompose_mass <- function(mass,
                                     ppm = 5,
                                     mzabs = NULL,
                                     elements = c("C", "H", "N", "O", "P", "S"),
                                     min_elements = NULL,
                                     max_elements = NULL,
                                     check_rule = FALSE) {

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
      nm <- names(x)
      v <- as.integer(x)
      names(v) <- nm
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

  if (isTRUE(check_rule) && any(min_counts < 0L)) {
    message("`check_rule` disabled: negative min_elements are for replacement/difference formulas.")
    check_rule <- FALSE
  }

  # NULL == Inf: that side of the min() window is unused.
  as_window <- function(x, nm) {
    if (is.null(x)) return(Inf)
    x <- as.numeric(x)
    if (length(x) != 1L || is.na(x) || x < 0) {
      stop("`", nm, "` must be NULL, Inf, or a single non-negative number.")
    }
    x
  }
  ppm_win <- as_window(ppm, "ppm")
  mzabs_win <- as_window(mzabs, "mzabs")
  if (!is.finite(ppm_win) && !is.finite(mzabs_win)) {
    stop("At least one of `ppm` or `mzabs` must be a finite tolerance (`NULL`/`Inf` means unused).")
  }

  call_one <- function(m_target) {
    ppm_abs <- if (is.finite(ppm_win)) ppm_win * abs(m_target) * 1e-6 else Inf
    abs_error <- min(ppm_abs, mzabs_win)
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
#' @param mz Ion m/z (scalar or vector). If `charge = 0`, treated as neutral mass.
#' @param charge Integer charge `z` (e.g. `+1`, `+2`, `-1`). Use `0` for neutral input.
#' @param ppm Allowed deviation in ppm. `NULL` is the same as `Inf` (unused).
#'   Default `5`. Combined with `mzabs` by taking the **tighter** (smaller)
#'   window on the **neutral** mass used for MCP; they are not added. See
#'   **Mass tolerance**.
#' @param mzabs Allowed absolute deviation in Dalton. `NULL` is the same as
#'   `Inf` (unused). Default `NULL` (ppm-only). Combined with `ppm` by taking
#'   the tighter window (see **Mass tolerance**).
#' @param elements Character vector of allowed elements.
#' @param min_elements See [chemform_decompose_mass()].
#' @param max_elements See [chemform_decompose_mass()].
#' @param check_rule If `TRUE`, keep only candidates that pass
#'   [chemform_check_seven_golden_rules()] (passed through to
#'   [chemform_decompose_mass()]). Default `FALSE`. Disabled automatically
#'   when any `min_elements` count is negative (see [chemform_decompose_mass()]).
#'
#' @details
#' Neutral mass conversion when `charge != 0`:
#' `M = mz * abs(charge) + e * charge`, with `e = 0.00054857990943`.
#'
#' ## Mass tolerance
#' Same rule as [chemform_decompose_mass()]: `NULL` and `Inf` mean unused.
#' `ppm` is converted to Dalton in R, then the **minimum** of the two windows
#' is passed to C++:
#'
#' `abs_error = min(ppm * |M| * 1e-6, mzabs)`
#'
#' on the converted **neutral** mass `M`. Defaults are `ppm = 5` and
#' `mzabs = NULL` (5 ppm). Set `ppm = NULL` (or `Inf`) for an absolute window
#' only.
#'
#' @return A data.frame with columns `formula`, `exactmass`, `mz`, `ppm`, `charge`,
#'   and `mz_target`. Rows are sorted by increasing `abs(ppm)` vs the input m/z.
#' @export
#' @seealso [chemform_decompose_mass()], [chemform_mz()], [chemform_check_seven_golden_rules()]
chemform_decompose_mz <- function(mz,
                                   charge = 0,
                                   ppm = 5,
                                   mzabs = NULL,
                                   elements = c("C", "H", "N", "O", "P", "S"),
                                   min_elements = NULL,
                                   max_elements = NULL,
                                   check_rule = FALSE) {

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
      max_elements = max_elements,
      check_rule = check_rule
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
