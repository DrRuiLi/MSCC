formula_calculate_lc8 <- function(Formula1 = "C2H4O1S2P1",Formula2 = "N1H1O-1",sign = 1,Valid_formula = FALSE ){

  to.return <- lc8::my_calculate_formula(Formula1 ,Formula2 ,sign ,Valid_formula )
  ### when formula calculate result to 0, such as "CH - CH = NULL", it will return "NANA"
  to.return[to.return == "NANA"] <- "C0"
  return(to.return)


}




#' @title isotopes_pattern_enviPat
#' @description update of enviPat::isopattern, which do not output the formula of isotopologues
#' this function return isotopologues with abundance more than 0.01%
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
isotopes_pattern_enviPat <- function(chemform) {
  #chemform <- "C80[13]C3H33[2]H12"
  data("isotopes",package = "enviPat")
  data("elem_table",package = "lc8")
  isopat <-
    enviPat::isopattern(isotopes = isotopes,
                        chemforms = chemform,
                        verbose =F,
                        threshold = 0.01)[[1]]
  iso.matrix <- isopat[, 3:ncol(isopat)]

  ele <- colnames(iso.matrix)
  #ele <- str_replace(ele,"[0-9]+",paste0("[",str_match(ele,"[0-9]+"),"]"))

  formulat_list <- c()
  for (i in 1:nrow(iso.matrix)) {
    formulat_list[i] <-
      data.frame(element = ele ,
                 n = iso.matrix[i, ]) %>%
      mutate(element = elem_table$element[match(element , elem_table$isotope)]) %>%
      group_by(element) %>%
      summarise(sum(n)) %>%
      filter(`sum(n)`!=0)%>%
      rowwise() %>%
      mutate(f = paste0(element, `sum(n)`)) %>%
      pull(f) %>%
      paste0(collapse = "")

  }
  formula_raw <- chemform
  select_elemet <- function(x ){
    ele_table <- lc8::my_break_formula(x)%>%as.data.frame()%>%
      filter(count >0)%>%
      mutate(f = paste0(elem,count))%>%
      pull(f)%>%
      paste0(collapse = "")
    ele_table
  }
  isopata <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])%>%
    rowwise()%>%
    mutate(isotope_element = formula_calculate_lc8(formula,formula_raw , -1),
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
chemform_mz_lc8 <- function(chemform = "C2H4O1S2P1"){

  lc8::formula_mz(chemform)

}



#' @title chemform_formate
#' @description
#' formate a chemical formula:
#' 1. replace `D` with `[1]H`
#'
#' @param chemform
#'
#' @return
#' @export
#'
#' @examples
chemform_formate <- function(chemform = "C11H22NO4"){

  enviPat::check_chemform(chemforms = chemform,isotopes=iso.table )$new_formula%>%
    sub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)### Replace D with [2]H


}




