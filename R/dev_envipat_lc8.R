#' @title chemical formula calculation
#' @description
#' two vector of `chemform`, return a matrix
#'
#'
#' @param Formula1
#' @param Formula2
#' @param sign
#' @param Valid_formula
#'
#' @return
#' @export
#'
#' @examples
chemform_calculate_matrix <- function(Formula1 = "C2H4O1S2P1",Formula2 = "N1H1O-1",sign = 1,Valid_formula = FALSE ){

  to.return <- lc8::my_calculate_formula(Formula1 ,Formula2 ,sign ,Valid_formula )
  ### when formula calculate result to 0, such as "CH - CH = NULL", it will return "NANA"
  to.return[to.return == "NANA"] <- NA
  return(to.return)


}

#' @title chemical formula calculation
#' @description
#' two vector of `chemform`, must be same length,return a vector
#'
#' @param Formula1
#' @param Formula2
#' @param sign
#'
#' @return
#' @export
#'
#' @examples
chemform_calculate_vector <- function(Formula1 = chem_formula_template  ,
                                      Formula2 = chem_formula_template ,
                                      sign = sample(1,replace = T,length(Formula1)) ){
  if (length(sign)==1 ) {
    sign <- rep(sign,length(Formula1))
  }
  any.na <- is.na(Formula1)|is.na(Formula2)|is.na(sign)

  to.return <- Formula1
  if (any(any.na)) {
    to.return[any.na] <-NA
   to.return[!any.na] <-chemform_calculate_vector(
    Formula1[!any.na],
    Formula2[!any.na] ,
    sign[!any.na]
    )

  }else{
    to.return <- sapply(1:length(Formula1),function(x){
      lc8::my_calculate_formula(Formula1[x] ,Formula2[x] ,sign[x] ,Valid_formula = F )
  })
  }

  to.return
  ### when formula calculate result to 0, such as "CH - CH = NULL", it will return "NANA"
  to.return[to.return == "NANA"] <- NA
  return(to.return)

}







#' @title chemform_isotopes_pattern_enviPat
#' @description update of enviPat::isopattern, which do not output the formula of isotopologues
#' this function return isotopologues with abundance more than 0.01%
#' @param chemform chemical formula
#'
#' @return
#' @export
#'
#' @examples
chemform_isotopes_pattern_enviPat <- function(chemform,thresh = 0.1) {
  #chemform <- "C80[13]C3H33[2]H12"
  #data("isotopes",package = "enviPat")
 # data("elem_table")
  elem_table <- MSCC::elem_table
  isopat <-
    enviPat::isopattern(isotopes = MSCC::isotopes,
                        chemforms = chemform_formate(chemform),
                        verbose =F,emass = 0.00054858,
                        threshold = thresh)[[1]]
  iso.matrix <- isopat[, 3:ncol(isopat)]
  if (is.null(nrow(iso.matrix))) {
    isopata <- tibble(formula = chemform,
                          m.z = chemform_mz_lc8(chemform),
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
      dplyr::mutate(element = elem_table$element[match(element , elem_table$isotope)]) %>%
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
    ele_table <- lc8::my_break_formula(x)%>%
      as.data.frame()%>%
      dplyr::filter(count >0)%>%
      dplyr::mutate(f = paste0(elem,count))%>%
      dplyr::pull(f)%>%
      paste0(collapse = "")
    ele_table
  }
  isopata <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])%>%
    dplyr::rowwise()%>%
    dplyr::mutate(isotope_element = chemform_calculate_matrix(formula,formula_raw , -1),
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
#' @return
#' @export
#'
#' @examples
chemform_mz_lc8 <- function(chemform = "C2H4O1S2P1",charge = 0){


  mz <- lc8::formula_mz( chemform,charge  = 0)### This function does not support charge as vector
  e_mass = 0.00054857990943
  mz= mz - e_mass * charge
  return(mz)
}




