chemform_calculate_lc8 <- function(Formula1 = "C2H4O1S2P1",Formula2 = "N1H1O-1",sign = 1,Valid_formula = FALSE ){

  to.return <- lc8::my_calculate_formula(Formula1 ,Formula2 ,sign ,Valid_formula )
  ### when formula calculate result to 0, such as "CH - CH = NULL", it will return "NANA"
  to.return[to.return == "NANA"] <- "C0"
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
chemform_isotopes_pattern_enviPat <- function(chemform) {
  #chemform <- "C80[13]C3H33[2]H12"
  #data("isotopes",package = "enviPat")
  #data("elem_table",package = "lc8")
  isopat <-
    enviPat::isopattern(isotopes = isotopes,
                        chemforms = chemform,
                        verbose =F,emass = 0.00054858,
                        threshold = 0.1)[[1]]
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
      mutate(element = elem_table$element[match(element , elem_table$isotope)]) %>%
      #mutate(element = paste0("[",element,"]")) %>%
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
    mutate(isotope_element = chemform_calculate_lc8(formula,formula_raw , -1),
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
#' 1. add 1 after single element
#' 2. replace `D` with `[1]H`
#'
#' @param chemform
#'
#' @return
#' @export
#'
#' @examples
chemform_formate <- function(chemform = "C11H22NO4D"){

  ### add 1 after elements
  chemform <- gsub(pattern = "([A-z](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
    gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)

  ### Replace D with [2]H
  #chemform <- enviPat::check_chemform(chemforms = chemform,isotopes=isotopes )$new_formula


  chemform
}



#' @title chemform_adduct
#' @description
#' get chemical formula with adduct
#'
#' @param chemform chemical formula
#' @param adduct adduct form, such as "[M+H]+"
#'
#' @return
#' @export
#'
#' @examples
chemform_adduct <- function(chemform = chem_formula_template,
                            adduct = "[M+H]+"){
  if ( !adduct%in% adduct.table$Adduct) {

    stop(paste0("Adduct form ",crayon::red(adduct), " incorrect, please check"))

  }
  chemform.diff <- adduct.table%>%
    dplyr::filter(Adduct == adduct)%>%
    dplyr::pull(Formula_diff)
  chemform_calculate_lc8(Formula1 = chemform,
                         Formula2 = chemform.diff,
                         sign = 1)%>%
    as.vector()




}


get.adduct.table.from.enviPat <- function(){


  data("adducts",package = "enviPat")

  adduct.table <- adducts%>%
    dplyr::rowwise()%>%
    dplyr::mutate(a=case_when(Formula_add=="FALSE"~"C0",
                              T~ Formula_add),
                  b=case_when(Formula_ded=="FALSE"~"C0",
                              T~ Formula_ded),
                  Formula_diff = chemform_calculate_lc8(
                    chemform_formate(a),
                    chemform_formate(b),
                    -1
                  ),
                  Adduct = paste0("[",Name,"]",ifelse(Ion_mode=="positive","+","-"))
    )%>%
    dplyr::select(-a,-b)
  return(adduct.table)
}
