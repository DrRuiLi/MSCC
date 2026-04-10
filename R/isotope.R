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
#' generate all possible isotope substitution combinations and their mass differences
#' for given element counts. Acceptable elements: C, H, O, N, P, S
#'
#' @param ... named arguments of element counts, e.g. `C = 3, N = 2, O = 5`
#'
#' @return data.table with columns `chemform_diff` and `mass_diff`
#' @export
#'
#' @examples get_isotope_mass_diff(C = 3, N = 2, O = 1)
get_isotope_mass_diff <- function(...) {
  args <- list(...)
  allowed <- c("C", "H", "O", "N", "P", "S")

  elements <- names(args)
  if (is.null(elements) || any(!elements %in% allowed)) {
    stop("All inputs must be named with allowed elements: ", paste(allowed, collapse = ", "))
  }

  counts <- as.integer(unlist(args))

  iso_map <- c(C = "[13]C", H = "[2]H", O = "[18]O", N = "[15]N", S = "[34]S")
  iso_minor <- iso_map[elements]
  has_iso <- !is.na(iso_minor)
  elements <- elements[has_iso]
  counts <- counts[has_iso]
  if (length(elements) == 0) {
    return(data.table::data.table(chemform_diff = character(0), mass_diff = numeric(0)))
  }
  iso_minor <- iso_minor[has_iso]
  mass_single <- mapply(isotope_mass_diff, iso_minor)

  ranges <- lapply(counts, function(n) seq_len(n + 1) - 1L)
  names(ranges) <- elements
  grid <- expand.grid(ranges)

  chemform_diff <- apply(grid, 1, function(row) {
    paste0(iso_minor, row, collapse = "")
  })

  mass_diff <- as.numeric(as.matrix(grid) %*% mass_single)

  data.table::data.table(chemform_diff = chemform_diff, mass_diff = mass_diff)
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
chemform_isotope_label <- function(chemform, ele, count) {
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
