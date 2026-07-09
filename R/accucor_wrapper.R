#' AccuCor-style natural isotope correction for ratio matrix
#'
#' @param raw_ratio Numeric matrix (rows = isotopologues, cols = samples).
#' @param formula Molecular formula string.
#' @param iso_ele Tracer isotope, e.g. `"[13]C"`, `"[2]H"`.
#' @param purity Tracer purity.
#' @param Resolution MS resolving power.
#' @param ResDefAt m/z where resolution is defined.
#'
#' @return Numeric matrix with corrected ratios.
#' @export
accucor_natural_correction <- function(raw_ratio,
                                       formula,
                                       iso_ele = "[13]C",
                                       purity = 1,
                                       Resolution = 140000,
                                       ResDefAt = 200) {
  if (is.null(raw_ratio) || !is.matrix(raw_ratio)) {
    stop("`raw_ratio` must be a numeric matrix.")
  }
  if (!is.numeric(raw_ratio)) {
    suppressWarnings(storage.mode(raw_ratio) <- "numeric")
  }
  if (!is.numeric(raw_ratio)) {
    stop("`raw_ratio` must be a numeric matrix.")
  }
  if (!nrow(raw_ratio) || !ncol(raw_ratio)) {
    return(raw_ratio)
  }
  if (missing(formula) || is.null(formula) || !nzchar(as.character(formula)[1])) {
    stop("`formula` must be provided.")
  }
  formula <- as.character(formula)[1]

  purity <- as.numeric(purity)[1]
  if (!is.finite(purity) || is.na(purity) || purity < 0 || purity > 1) {
    stop("`purity` must be a numeric value between 0 and 1.")
  }
  Resolution <- as.numeric(Resolution)[1]
  ResDefAt <- as.numeric(ResDefAt)[1]
  if (!is.finite(Resolution) || is.na(Resolution) || Resolution <= 0) {
    stop("`Resolution` must be a positive numeric value.")
  }
  if (!is.finite(ResDefAt) || is.na(ResDefAt) || ResDefAt <= 0) {
    stop("`ResDefAt` must be a positive numeric value.")
  }

  iso_ele <- as.character(iso_ele)[1]
  if (identical(iso_ele, "D")) iso_ele <- "[2]H"
  label_idx <- suppressWarnings(as.integer(gsub(".*?([0-9]+)$", "\\1", rownames(raw_ratio))))
  if (all(is.na(label_idx))) label_idx <- seq_len(nrow(raw_ratio)) - 1L
  label_idx[is.na(label_idx)] <- seq_len(sum(is.na(label_idx))) - 1L

  atom_counts <- tryCatch(CHNOSZ::makeup(formula), error = function(e) NULL)
  if (is.null(atom_counts)) {
    stop("Unable to parse `formula`.")
  }

  if (identical(iso_ele, "[13]C")) {
    corrected_full <- accucor::carbon_isotope_correction(
      formula = formula,
      datamatrix = raw_ratio,
      label = label_idx,
      Resolution = Resolution,
      ResDefAt = ResDefAt,
      purity = purity
    )
    corrected <- corrected_full[label_idx + 1L, , drop = FALSE]
    dimnames(corrected) <- dimnames(raw_ratio)
    corrected[is.na(corrected)] <- 0
    corrected[corrected < 0] <- 0
    return(corrected)
  }

  elem_table <- get_elem_table()
  if (!is.data.frame(elem_table) ||
      !all(c("element", "abundance") %in% colnames(elem_table))) {
    stop("`MSCC::elem_table` does not contain expected columns: element, abundance.")
  }
  tracer_row <- elem_table[as.character(elem_table$element) %in% iso_ele, , drop = FALSE]
  if (!nrow(tracer_row)) {
    stop("`iso_ele` not found in `MSCC::elem_table$element`.")
  }
  r <- as.numeric(tracer_row$abundance[[1]])
  if (!is.finite(r) || is.na(r) || r < 0) {
    stop("Invalid abundance for `iso_ele` in `MSCC::elem_table`.")
  }
  natural_p <- r / (1 + r)

  tracer_element <- sub("^\\[[0-9]+\\]", "", iso_ele)
  if (!nzchar(tracer_element)) tracer_element <- iso_ele
  if (is.na(atom_counts[tracer_element])) {
    stop("Tracer element is missing in `formula`.")
  }
  tracer_mass_diff <- as.numeric(tracer_row$Mass_Dif[[1]])
  if (!is.finite(tracer_mass_diff) || is.na(tracer_mass_diff) || tracer_mass_diff <= 0) {
    stop("Invalid tracer `Mass_Dif` in `MSCC::elem_table`.")
  }

  n_tracer <- as.integer(atom_counts[tracer_element])
  if (max(label_idx, na.rm = TRUE) > n_tracer) {
    stop("Observed isotopologue rows exceed tracer atom count from `formula`.")
  }
  get_atom_n <- function(ele) {
    v <- atom_counts[ele]
    if (length(v) == 0 || is.na(v)) return(0L)
    as.integer(v[[1]])
  }
  labels <- 0:n_tracer
  exp_matrix <- matrix(0, nrow = n_tracer + 1L, ncol = ncol(raw_ratio))
  for (i in seq_len(nrow(raw_ratio))) {
    idx <- label_idx[i] + 1L
    if (idx >= 1L && idx <= nrow(exp_matrix)) exp_matrix[idx, ] <- as.numeric(raw_ratio[i, ])
  }
  exp_matrix[is.na(exp_matrix)] <- 0

  purity_matrix <- diag(n_tracer + 1L)
  for (i in seq_len(n_tracer + 1L)) {
    purity_matrix[i, ] <- vapply(labels, function(x) {
      stats::dbinom(x - i + 1L, x, (1 - purity))
    }, numeric(1))
  }

  tracer_matrix <- matrix(0, nrow = n_tracer + 1L, ncol = n_tracer + 1L)
  for (i in seq_len(n_tracer + 1L)) {
    tracer_matrix[, i] <- vapply(labels, function(x) {
      stats::dbinom(x - i + 1L, n_tracer - i + 1L, natural_p)
    }, numeric(1))
  }

  iso_defs <- data.frame(
    iso = c("[13]C", "[2]H", "[15]N", "[17]O", "[18]O", "[33]S", "[34]S", "[29]Si", "[30]Si", "[37]Cl", "[81]Br"),
    base = c("C", "H", "N", "O", "O", "S", "S", "Si", "Si", "Cl", "Br"),
    weight = c(1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  iso_defs <- iso_defs[iso_defs$iso != iso_ele, , drop = FALSE]
  iso_defs$atom_n <- vapply(iso_defs$base, get_atom_n, integer(1))

  get_iso_prob <- function(iso_key, base_ele) {
    x <- elem_table[as.character(elem_table$element) %in% c(base_ele, iso_key), , drop = FALSE]
    if (!nrow(x)) return(NA_real_)
    el <- as.character(x$element)
    ab <- as.numeric(x$abundance)
    keep <- el %in% c(base_ele, iso_key) & is.finite(ab)
    el <- el[keep]
    ab <- ab[keep]
    if (!length(ab) || !any(el %in% iso_key)) return(NA_real_)
    ab[match(iso_key, el)] / sum(ab)
  }
  get_iso_massdiff <- function(iso_key) {
    x <- elem_table[as.character(elem_table$element) %in% iso_key, , drop = FALSE]
    if (!nrow(x) || !("Mass_Dif" %in% colnames(x))) return(NA_real_)
    as.numeric(x$Mass_Dif[[1]])
  }

  iso_defs$prob <- mapply(get_iso_prob, iso_defs$iso, iso_defs$base, SIMPLIFY = TRUE)
  iso_defs$mass_diff <- vapply(iso_defs$iso, get_iso_massdiff, numeric(1))
  iso_defs <- iso_defs[
    iso_defs$atom_n > 0 &
      is.finite(iso_defs$prob) &
      iso_defs$prob > 0 &
      is.finite(iso_defs$mass_diff),
    ,
    drop = FALSE
  ]

  get_mono_mass <- function(ele) {
    x <- elem_table[as.character(elem_table$element) %in% ele, , drop = FALSE]
    if (!nrow(x)) return(NA_real_)
    as.numeric(x$mass[[1]])
  }
  molecular_weight <- sum(vapply(names(atom_counts), function(ele) {
    mm <- get_mono_mass(ele)
    if (!is.finite(mm) || is.na(mm)) return(0)
    mm * as.numeric(atom_counts[[ele]])
  }, numeric(1)), na.rm = TRUE)
  mass_limit <- 1.66 * molecular_weight^1.5 / Resolution / sqrt(ResDefAt)

  nontracer_matrix <- matrix(0, nrow = n_tracer + 1L, ncol = n_tracer + 1L)
  if (!nrow(iso_defs)) {
    diag(nontracer_matrix) <- 1
  } else {
    grid_list <- lapply(seq_len(nrow(iso_defs)), function(i) 0:iso_defs$atom_n[i])
    names(grid_list) <- iso_defs$iso
    combos <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    if (nrow(combos)) {
      by_base <- split(iso_defs$iso, iso_defs$base)
      for (base_ele in names(by_base)) {
        keys <- by_base[[base_ele]]
        if (length(keys) > 1L) {
          combos <- combos[rowSums(combos[, keys, drop = FALSE]) <= get_atom_n(base_ele), , drop = FALSE]
        }
      }
    }

    if (nrow(combos)) {
      mass_sum <- rep(0, nrow(combos))
      mass_delta <- rep(0, nrow(combos))
      for (k in seq_len(nrow(iso_defs))) {
        key <- iso_defs$iso[k]
        cnt <- as.numeric(combos[[key]])
        mass_sum <- mass_sum + cnt * iso_defs$weight[k]
        mass_delta <- mass_delta + cnt * iso_defs$mass_diff[k]
      }
      keep <- mass_sum <= n_tracer & abs(mass_delta - tracer_mass_diff * mass_sum) < mass_limit
      combos <- combos[keep, , drop = FALSE]
      mass_sum <- mass_sum[keep]

      calc_prob <- function(one_row) {
        p <- 1
        for (base_ele in names(by_base)) {
          keys <- by_base[[base_ele]]
          n_base <- get_atom_n(base_ele)
          counts <- as.numeric(one_row[keys])
          probs <- iso_defs$prob[match(keys, iso_defs$iso)]
          if (length(keys) == 1L) {
            p <- p * stats::dbinom(counts[[1]], n_base, probs[[1]])
          } else {
            p0 <- 1 - sum(probs)
            if (p0 < 0) p0 <- 0
            p <- p * stats::dmultinom(c(n_base - sum(counts), counts), n_base, c(p0, probs))
          }
        }
        p
      }
      probs <- apply(combos, 1, calc_prob)
      for (ii in seq_len(nrow(combos))) {
        ms <- as.integer(mass_sum[ii])
        if (ms < 0 || ms > n_tracer) next
        for (jj in seq_len(n_tracer + 1L - ms)) {
          nontracer_matrix[ms + jj, jj] <- nontracer_matrix[ms + jj, jj] + probs[ii]
        }
      }
    }
    if (!any(nontracer_matrix > 0)) diag(nontracer_matrix) <- 1
  }

  corrected_full <- matrix(0, nrow = n_tracer + 1L, ncol = ncol(raw_ratio))
  for (j in seq_len(ncol(exp_matrix))) {
    y <- as.numeric(exp_matrix[, j])
    y[is.na(y)] <- 0
    fit <- nnls::nnls(nontracer_matrix %*% tracer_matrix %*% purity_matrix, y)
    corrected_full[, j] <- stats::coef(fit)
  }

  corrected <- corrected_full[label_idx + 1L, , drop = FALSE]
  dimnames(corrected) <- dimnames(raw_ratio)
  corrected[is.na(corrected)] <- 0
  corrected[corrected < 0] <- 0
  corrected
}
