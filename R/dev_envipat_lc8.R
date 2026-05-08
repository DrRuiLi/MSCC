
#' @title chemform_isotopes_pattern_enviPat
#' @description update of enviPat::isopattern, which do not output the formula of isotopologues
#' this function return isotopologues with abundance more than 0.01%
#' @param chemform chemical formula
#'
#' @return isotopes pattern
#' @export
#'

chemform_isotopes_pattern_enviPat <- function(chemform,thresh = 0.1) {
  #chemform <- "C80[13]C3H33[2]H12"
  #data("isotopes",package = "enviPat")
 # data("elem_table")
  elem_table <- MSCC::elem_table
  isopat <-
    enviPat::isopattern(isotopes = MSCC::isotopes_from_envipat,
                        chemforms = chemform_formate(chemform),
                        verbose =F,emass = 0.00054857990924,
                        threshold = thresh)[[1]]
  iso.matrix <- isopat[, 3:ncol(isopat)]
  if (is.null(nrow(iso.matrix))) {
    isopata <- tibble::tibble(formula = chemform,
                          m.z = chemform_mz(chemform),
                          abundance =100,
                          isotope_element = ""

                          )
    return(isopata)

  }

  ele <- colnames(iso.matrix)
  #ele <- str_replace(ele,"[0-9]+",paste0("[",str_match(ele,"[0-9]+"),"]"))

  formulat_list <- c()
  for (i in 1:nrow(iso.matrix)) {
    formulat_list[i] <-
      data.frame(element = ele ,
                 n = iso.matrix[i, ]) %>%
      dplyr::mutate(element = MSCC::elem_table$element[match(element , MSCC::elem_table$isotope)]) %>%
      #mutate(element = paste0("[",element,"]")) %>%
      dplyr::group_by(element) %>%
      dplyr::summarise(sum(n)) %>%
      dplyr::filter(`sum(n)`!=0)%>%
      dplyr::rowwise() %>%
      dplyr::mutate(f = paste0(element, `sum(n)`)) %>%
      dplyr::pull(f) %>%
      paste0(collapse = "")

  }
  formula_raw <- chemform
  select_elemet <- function(x ){
    if (is.na(x)) {
      return(NA)
    }
    x.ele.m <- chemform_parse(x)
    x.ele.m[x.ele.m<0] <- 0
    chemform_from_ele_matrix(x.ele.m)
  }
  isopata <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])%>%
    dplyr::rowwise()%>%
    dplyr::mutate(
      isotope_element = chemform_calc(formula,formula_raw , "-",return = "chemform"),
      isotope_element = select_elemet(isotope_element))%>%
    dplyr::arrange(-abundance)


  return(isopata)

}





#' @title chemform_mz_lc8
#' @description
#' calculate mz of a given chemical formula
#'
#' @param chemform, such as `"C2H4"`
#'
#' @return mz
#' @export
#'

chemform_mz_lc8 <- function(chemform = "C2H4O1S2P1",charge = 0){


  mz <- lc8::formula_mz( chemform,charge  = 0)### This function does not support charge as vector
  e_mass = 0.00054857990943
  mz= mz - e_mass * charge
  return(mz)
}




