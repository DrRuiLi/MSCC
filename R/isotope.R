#' @title isotope_mass_diff
#' @description
#' calculate the mass difference of a given isotope to its major (most abundant) form
#'
#' @param isotope string, isotope notation such as `\\[13\\]C`, `\\[2\\]H`
#'
#' @return numeric, mass difference (M\\[13\\]C - M\\[12\\]C)
#' @export
#'
#' @examples isotope_mass_diff("[13]C")
#' isotope_mass_diff("[2]H")
isotope_mass_diff <- function(isotope) {
  iso_data <- MSCC::elem_table

  mass_num <- stringr::str_extract(isotope, "[0-9]+")
  element <- stringr::str_extract(isotope, "[A-Z][a-z]?")
  iso_name <- paste0(mass_num, element)

  iso_mass <- iso_data$mass[match(iso_name, iso_data$isotope)]
  if (is.na(iso_mass)) {
    stop(paste0("Isotope not found: ", isotope))
  }

  ele_rows <- iso_data[iso_data$element == element, ]
  major_mass <- ele_rows$mass[which.max(ele_rows$abundance)]

  return(iso_mass - major_mass)
}


#' @title get_isotope_mass_diff
#' @description
#' For each formula / element fragment in \code{element}, runs
#' \code{chemform_isotopes_pattern_enviPat()} and returns a named numeric
#' vector of isotopologue mass differences vs the monoisotopic peak.
#' Names use compact isotope notation (e.g. \[13\]C, \[13\]C2, \[34\]S).
#' Accepts plain element symbols (\code{"C"}, \code{"S"}) or counted
#' fragments (\code{"C10"}).
#'
#' @param element Character vector of formulas / element symbols,
#'   e.g. \code{c("C10", "H10", "O2", "S")}.
#' @param threshold enviPat abundance cutoff passed as \code{thresh}
#'   (percent of the monoisotopic peak; default \code{0.0001}).
#'
#' @return Named numeric vector: names are isotope labels, values are
#'   \code{mass_diff}.
#' @export
#'
#' @examples
#' get_isotope_mass_diff(c("C10", "H10", "O2", "S"))
#' get_isotope_mass_diff(element = c("C", "N", "S"), threshold = 0.01)
get_isotope_mass_diff <- function(element = c("C10","H10","O5","N5","P3","S3","K","Cl","Br"), threshold = 0.0001) {
  if ( is.null(element)) {
    stop("`element` must be a character vector of formulas / element symbols")
  }
  if (!is.character(element)) {
    stop("`element` must be a character vector of formulas / element symbols")
  }
  threshold <- as.numeric(threshold)[[1]]
  if (!is.finite(threshold) || threshold < 0) {
    stop("`threshold` must be a non-negative number")
  }
  .isotope_mass_diff_for_elements(element, threshold = threshold)
}

#' Compact enviPat isotope_element labels: `[13]C1` -> `[13]C`, keep `[13]C2`.
#' @keywords internal
.isotope_label_compact <- function(x) {
  gsub("(\\[[0-9]+\\][A-Za-z]+)1(?![0-9])", "\\1", as.character(x), perl = TRUE)
}

#' @keywords internal
.isotope_mass_diff_for_elements <- function(elements, threshold = 0.0001) {
  empty <- setNames(numeric(0), character(0))
  elements <- as.character(elements)
  elements <- elements[nzchar(elements) & !is.na(elements)]
  if (!length(elements)) {
    return(empty)
  }

  out_names <- character(0)
  out_vals <- numeric(0)

  for (fm in elements) {
    pat <- tryCatch(
      chemform_isotopes_pattern_enviPat(fm, thresh = threshold),
      error = function(e) {
        warning("enviPat pattern failed for '", fm, "': ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(pat) || !nrow(pat)) {
      next
    }
    pat <- as.data.frame(pat, stringsAsFactors = FALSE)
    iso <- as.character(pat$isotope_element)
    mono_idx <- which(is.na(iso) | !nzchar(iso))
    if (!length(mono_idx)) {
      mono_mz <- chemform_mz(fm)
    } else {
      mono_mz <- as.numeric(pat$m.z[[mono_idx[[1]]]])
    }
    keep <- !(is.na(iso) | !nzchar(iso))
    if (!any(keep)) {
      next
    }
    labs <- .isotope_label_compact(iso[keep])
    vals <- as.numeric(pat$m.z[keep]) - mono_mz
    # keep first occurrence if a label repeats across fragments
    new <- !labs %in% out_names
    out_names <- c(out_names, labs[new])
    out_vals <- c(out_vals, vals[new])
  }

  setNames(out_vals, out_names)
}


#' @title chemform_isotope_label
#' @description
#' replace a specified number of atoms of an element with an isotope in a chemical formula
#'
#' @param chemform chemical formula
#' @param ele isotope notation, e.g. `\\[13\\]C`, `\\[2\\]H`
#' @param count number of atoms to replace with the isotope
#'
#' @return character, the labeled chemical formula
#' @export
#'
#' @examples chemform_isotope_label("C6H12O6", "[13]C", 3)
chemform_isotope_label <- function(chemform = demo_chemform,
                                    ele = "[13]C",
                                    count = 1) {
  mass_num <- stringr::str_extract(ele, "[0-9]+")
  element  <- stringr::str_extract(ele, "[A-Z][a-z]?")

  if (is.na(mass_num) || is.na(element)) {
    stop("Invalid isotope format: ", ele)
  }

  chemform.formated <- chemform_formate(chemform)
  chemform.matrix <- chemform_parse(chemform.formated, return = "matrix")

  iso_map <- c(C = "[13]C", H = "[2]H", O = "[18]O", N = "[15]N", S = "[34]S")
  iso_minor <- iso_map[element]
  if (is.na(iso_minor)) {
    iso_minor <- ele
  }

  sapply(seq_len(nrow(chemform.matrix)), function(i) {
    ele_count <- chemform.matrix[i, element, drop = TRUE]
    if (is.na(ele_count)) ele_count <- 0
    actual_count <- min(count, ele_count)
    chemform.matrix[i, element] <- ele_count - actual_count
    chemform.major <- chemform_from_ele_matrix(chemform.matrix[i, , drop = FALSE])[1]
    chemform.matrix[i, element] <- ele_count
    if (actual_count > 0) {
      paste0(iso_minor, actual_count, chemform.major)
    } else {
      chemform.major
    }
  })
}
