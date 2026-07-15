# Internal schema harmonization for legacy MSdev chemistry callers.

#' Zero-pad integers to fixed-width strings (MSCC-local; do not rely on MSdev)
#' @keywords internal
num2str <- function(x, n.digit = NA) {
  if (length(x) == 0L) return(character(0))
  if (is.na(n.digit)) n.digit <- max(nchar(as.character(x)), na.rm = TRUE)
  sprintf(paste0("%0", n.digit, "d"), as.integer(x))
}

#' Extract first numeric token from strings (MSCC-local; do not rely on MSdev)
#' @keywords internal
str_extract_num <- function(x) {
  as.numeric(stringr::str_extract(x, "[:digit:]+"))
}

#' Named vector helper (MSCC-local; do not rely on MSdev)
#' @keywords internal
make_vector <- function(x = NA, name = NULL) {
  if (length(x) == length(name)) {
    names(x) <- name
  }
  if (length(x) == 1L && length(name)) {
    x <- rep(x, length(name))
    names(x) <- name
  }
  x
}

#' Subset igraph to selected vertices (MSCC-local)
#' @keywords internal
igraph_filter_vertex <- function(ig, v) {
  if (is.numeric(v) || is.logical(v) || is.character(v)) {
    v <- igraph::V(ig)[v]
  }
  igraph::delete.vertices(ig, setdiff(names(igraph::V(ig)), names(v)))
}

#' Keep vertices within graph distance of seeds (MSCC-local)
#' @keywords internal
igraph_filter_distance <- function(ig, from, dis = 1, ...) {
  dis.matrix <- igraph::distances(ig, from)
  dis.pass <- dis.matrix <= dis
  id <- apply(dis.pass, 2, function(x) any(x))
  igraph_filter_vertex(ig, id)
}

#' Index a matrix by rownames/colnames with NA fill (MSCC-local)
#' @keywords internal
get_matrix_value_fill_with_NA <- function(mat,
                                          rownames_vec = rownames(mat),
                                          colnames_vec = colnames(mat),
                                          drop = TRUE) {
  result_matrix <- matrix(NA,
                          nrow = length(rownames_vec),
                          ncol = length(colnames_vec))
  rownames(result_matrix) <- rownames_vec
  colnames(result_matrix) <- colnames_vec
  for (i in seq_along(rownames_vec)) {
    for (j in seq_along(colnames_vec)) {
      rowname <- rownames_vec[i]
      colname <- colnames_vec[j]
      if (rowname %in% rownames(mat) && colname %in% colnames(mat)) {
        result_matrix[i, j] <- mat[rowname, colname]
      }
    }
  }
  if (drop && length(result_matrix) == 1L) {
    result_matrix <- as.vector(result_matrix)
  }
  result_matrix
}

.ensure_elem_table_schema <- function() {
  et <- MSCC::elem_table
  if (!is.data.frame(et)) return(et)
  if (!"symbol" %in% colnames(et)) {
    et$symbol <- get_ele_uniso(as.character(et$element))
  }
  if (!"is.isotope" %in% colnames(et)) {
    et$is.isotope <- grepl("^\\[[0-9]+\\]", as.character(et$element))
  }
  if (!"Mass_Dif" %in% colnames(et)) {
    et$Mass_Dif <- 0
    iso.idx <- which(et$is.isotope)
    if (length(iso.idx)) {
      et$Mass_Dif[iso.idx] <- vapply(et$element[iso.idx], function(x) {
        tryCatch(isotope_mass_diff(x), error = function(e) NA_real_)
      }, numeric(1))
    }
  }
  et
}

#' Get harmonized element table
#'
#' @return Element metadata table with legacy compatibility columns.
#' @export
get_elem_table <- function() {
  .ensure_elem_table_schema()
}

#' Convert isotope notation to element symbol
#'
#' @param iso_ele Isotope notation such as `"[13]C"` or `"C"`.
#'
#' @return Element symbol without isotope mass prefix.
#' @export
get_ele_uniso <- function(iso_ele = "[13]C") {
  stringr::str_replace(as.character(iso_ele), "^\\[[0-9]+\\]", "")
}

#' Convert element to canonical isotope notation
#'
#' @param iso_ele Element symbol or isotope notation.
#'
#' @return Canonical isotope notation for supported tracer elements.
#' @export
trans_iso_ele <- function(iso_ele = "[13]C") {
  ele <- get_ele_uniso(iso_ele)
  switch(ele,
    C = "[13]C",
    H = "[2]H",
    N = "[15]N",
    O = "[18]O",
    S = "[34]S",
    P = "[31]P",
    iso_ele
  )
}

# Internal utility retained here to avoid circular dependency on MSdev.
list2df <- function(x) {
  if (is.null(x) || !length(x)) {
    return(data.frame())
  }
  x.names <- unique(unlist(lapply(x, names)))
  out <- matrix(NA_real_, nrow = length(x), ncol = length(x.names))
  colnames(out) <- x.names
  rownames(out) <- names(x)
  for (i in seq_along(x)) {
    xi <- x[[i]]
    if (is.null(xi) || !length(xi)) next
    idx <- match(names(xi), x.names)
    out[i, idx[!is.na(idx)]] <- as.numeric(xi[!is.na(idx)])
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Count an element in formula
#'
#' @param formula Chemical formula.
#' @param ele Element symbol or isotope notation.
#'
#' @return Integer vector of element counts.
#' @export
get_formula_ele_count <- function(formula, ele = "C") {
  ele <- get_ele_uniso(ele)
  mat <- chemform_parse(formula, return = "matrix")
  if (!ele %in% colnames(mat)) {
    return(rep(0L, nrow(mat)))
  }
  as.integer(mat[, ele, drop = TRUE])
}

#' Build isotopologue natural ratio baseline
#'
#' @param formula Chemical formula.
#' @param iso_ele Isotope notation.
#' @param ratio_matrix Numeric matrix with isotopologues as rows.
#'
#' @return Matrix with natural baseline ratios aligned to `ratio_matrix`.
#' @export
get_iso_natural_ratio <- function(formula, iso_ele, ratio_matrix) {
  if (is.null(ratio_matrix) || !is.matrix(ratio_matrix) || !nrow(ratio_matrix)) {
    return(ratio_matrix)
  }
  iso_table <- chemform_isotopes_pattern_enviPat(formula)
  iso_form <- as.character(iso_table$isotope_element)
  iso_count <- get_formula_ele_count(iso_form, ele = iso_ele)
  iso_prob <- as.numeric(iso_table$abundance) / 100
  natural <- setNames(iso_prob, iso_count)
  row_iso <- suppressWarnings(as.integer(gsub(".*?([0-9]+)$", "\\1", rownames(ratio_matrix))))
  if (all(is.na(row_iso))) {
    row_iso <- seq_len(nrow(ratio_matrix)) - 1L
  }
  out <- matrix(0, nrow = nrow(ratio_matrix), ncol = ncol(ratio_matrix), dimnames = dimnames(ratio_matrix))
  for (i in seq_len(nrow(out))) {
    k <- as.character(row_iso[i])
    if (!is.na(k) && k %in% names(natural)) {
      out[i, ] <- natural[[k]]
    }
  }
  out
}

#' Compute pairwise adduct mass differences by polarity
#'
#' @param pol Polarity (`0`/`1`, `"negative"`/`"positive"`).
#'
#' @return Data frame with adduct pairs and delta mass.
#' @export
get_adduct_mass_diff <- function(pol = 0) {
  pol_chr <- as.character(pol)
  ion_mode <- if (pol_chr %in% c("1", "positive", "Positive")) "positive" else "negative"
  adt <- MSCC::adduct.table
  adt <- adt[as.character(adt$Ion_mode) %in% ion_mode, , drop = FALSE]
  if (!nrow(adt)) return(data.frame())
  expand.grid(
    adduct1 = as.character(adt$Adduct),
    adduct2 = as.character(adt$Adduct),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      mass1 = adt$Mass[match(adduct1, adt$Adduct)],
      mass2 = adt$Mass[match(adduct2, adt$Adduct)],
      mass_diff = mass2 - mass1
    )
}

#' Infer polarity from adduct expression
#'
#' @param adduct Character adduct annotation, e.g. `"[M+H]+"`, `"[M-H]-"`.
#'
#' @return Integer polarity (`0` negative, `1` positive).
#' @export
get_polarity_from_adduct <- function(adduct) {
  if (is.numeric(adduct) || is.integer(adduct)) {
    if (adduct == 0) return(0L)
    if (adduct == 1) return(1L)
  }
  ad <- as.character(adduct)[[1]]
  if (identical(ad, "0") || identical(ad, "-")) return(0L)
  if (identical(ad, "1") || identical(ad, "+")) return(1L)
  if (grepl("\\[M\\-H\\]", ad) || grepl("\\[M\\-\\]", ad)) {
    return(0L)
  }
  1L
}

#' Match theoretical isotopes to xcms feature definitions
#'
#' @param isotopes_table Data frame from `chemform_isotopes_pattern_enviPat()`.
#' @param featuredef xcms feature definition data frame.
#' @param mz.ppm m/z tolerance in ppm.
#' @param rt.tol RT tolerance in seconds if both sides have RT.
#'
#' @return Matched isotope-feature table.
#' @export
match_isotopes_to_featuredef <- function(isotopes_table, featuredef, mz.ppm = 10, rt.tol = 10) {
  if (is.null(isotopes_table) || !nrow(isotopes_table) || is.null(featuredef) || !nrow(featuredef)) {
    return(data.frame())
  }
  iso <- as.data.frame(isotopes_table, stringsAsFactors = FALSE)
  fdf <- as.data.frame(featuredef, stringsAsFactors = FALSE)
  if (!"m.z" %in% colnames(iso) || !"mzmed" %in% colnames(fdf)) return(data.frame())
  if (!"feature_id" %in% colnames(fdf)) fdf$feature_id <- rownames(fdf)
  if (!"rtmed" %in% colnames(fdf)) fdf$rtmed <- NA_real_
  if (!"rt" %in% colnames(iso)) iso$rt <- NA_real_

  out <- lapply(seq_len(nrow(iso)), function(i) {
    mz <- suppressWarnings(as.numeric(iso$m.z[i]))
    rt <- suppressWarnings(as.numeric(iso$rt[i]))
    if (!is.finite(mz)) return(NULL)
    mz_err <- abs(suppressWarnings(as.numeric(fdf$mzmed)) - mz) / mz
    hit <- mz_err <= mz.ppm * 1e-6
    if (is.finite(rt)) {
      rt_err <- abs(suppressWarnings(as.numeric(fdf$rtmed)) - rt)
      hit <- hit & (rt_err <= rt.tol)
    }
    idx <- which(hit)
    if (!length(idx)) return(NULL)
    cbind(
      iso[rep(i, length(idx)), , drop = FALSE],
      fdf[idx, , drop = FALSE],
      mz_error_ppm = mz_err[idx] * 1e6,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame())
  dplyr::bind_rows(out)
}

#' Append feature intensity matrix to isotope matches
#'
#' @param matched_table Output from `match_isotopes_to_featuredef()`.
#' @param featureval Numeric matrix from `xcms::featureValues()`.
#'
#' @return Data frame with matched rows and joined intensity columns.
#' @export
match_isotopes_to_featureval <- function(matched_table, featureval) {
  if (is.null(matched_table) || !nrow(matched_table)) return(data.frame())
  out <- as.data.frame(matched_table, stringsAsFactors = FALSE)
  if (is.null(featureval) || !is.matrix(featureval) || !nrow(featureval)) return(out)
  fv <- as.data.frame(featureval, stringsAsFactors = FALSE)
  fv$feature_id <- rownames(featureval)
  val_cols <- setdiff(colnames(fv), "feature_id")
  colnames(fv)[match(val_cols, colnames(fv))] <- paste0("featureval_", val_cols)
  dplyr::left_join(out, fv, by = "feature_id")
}
