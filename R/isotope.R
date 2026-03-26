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
