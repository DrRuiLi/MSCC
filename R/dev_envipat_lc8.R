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
  data("elem_table")
  isopat <-
    enviPat::isopattern(isotopes = MSCC::isotopes,
                        chemforms = chemform,
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
    ele_table <- lc8::my_break_formula(x)%>%as.data.frame()%>%
      dplyr::filter(count >0)%>%
      dplyr::mutate(f = paste0(elem,count))%>%
      dplyr::pull(f)%>%
      paste0(collapse = "")
    ele_table
  }
  isopata <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])%>%
    dplyr::rowwise()%>%
    dplyr::mutate(isotope_element = chemform_calculate_lc8(formula,formula_raw , -1),
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
  chemform <- gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
    gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)


  ### remove error
  chemform[is.na(chemform)] <-"NA"###
  idx.error <- enviPat::check_chemform(isotopes = MSCC::isotopes,
                                       chemforms = chemform)
  chemform[which(idx.error$warning)] <-NA

  return(chemform)
  ### Replace D with [2]H
  #chemform <- enviPat::check_chemform(chemforms = chemform,isotopes=isotopes )$new_formula


  #chemform
}



#' @title chemform_adduct_formula
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

  #chemform <- chemform_formate(chemform)

  if (length(adduct)==1 ) {
    adduct <- rep(adduct,length(chemform))
  }
  #adduct <- sample(MSCC::adduct.table$Adduct,replace = T,length(chemform))

  adduct.check <- chemform_adduct_check(adduct)


  chem_df <- data.frame(chemform = chemform,
                        adduct = adduct.check$adduct.formated)%>%
    dplyr::mutate(chemform = chemform_formate(chemform),
                  MSCC::adduct.table[match(adduct,MSCC::adduct.table$Adduct),c("Formula_diff","Multi","Charge")],
                  chemform.adduct = chemform_calculate_vector(chemform,Formula_diff),
                  chemform.adduct.mass = chemform_mz_lc8(chemform.adduct,charge = Charge),
                  chemform.adduct.mz = chemform.adduct.mass*Multi/abs(Charge)
                  )





  return(chem_df)


}


get.adduct.table.from.enviPat <- function(){


  data("adducts",package = "enviPat")

  adduct.table <- adducts%>%
    dplyr::add_row(Name = "M-H2O+H",calc = "M-14.987089588",
                   Charge=-1,Mult=1,Mass=chemform_mz_lc8("H-1O-1"),Ion_mode="positive",
                   Formula_add="FALSE",Formula_ded="H1O1",Multi=1)%>%
    dplyr::filter(Name != "2M+3H2O+2H")%>%
    dplyr::rowwise()%>%
    dplyr::mutate(a=case_when(Formula_add=="FALSE"~"C0",
                              T~ Formula_add),
                  b=case_when(Formula_ded=="FALSE"~"C0",
                              T~ Formula_ded),
                  Formula_diff = chemform_calculate_vector(
                    chemform_formate(a),
                    chemform_formate(b),
                    -1
                  ),
                  Adduct = paste0("[",Name,"]",ifelse(Ion_mode=="positive","+","-"))
    )%>%
    dplyr::select(-a,-b)

  ### custom
  {
    adduct.table <- adduct.table %>%
      dplyr::mutate(Adduct_Syn = paste0(Name,";",Adduct,";"),
                    Adduct_Syn = case_when(
                    Name == "M+Hac-H" ~paste0( Adduct_Syn,"[M+CH3COO]-;[M+CH3COOH-H]-;"),
                    Name == "2M+Hac-H" ~paste0( Adduct_Syn,"[2M+CH3COO]-;"),
                    Name == "M+FA-H" ~paste0( Adduct_Syn,"[M+HCOO]-;"),
                    Name == "M+" ~paste0( Adduct_Syn,"[M]+;"),
                    Name == "M-" ~paste0( Adduct_Syn,"[M]-;"),
                    Name == "M-H2O+H"~paste0( Adduct_Syn,"[M+H-H2O]+;"),


                      T~Adduct_Syn
                    ))

    }

 #  use_data(adduct.table,overwrite = T)

  return(adduct.table)
}



#' @chemform_adduct_check
#' @description
#' check if string is a adduct (match in adduct.table)
#'
#'
#' @param adduct.to.check
#'
#' @return
#' @export
#'
#' @examples
chemform_adduct_check <- function(adduct.to.check ){

  .check_adduct <- function(x){

    grepl(x = MSCC::adduct.table$Adduct_Syn,
          pattern = paste0(x,";"),
          fixed = T)->x.exist
    adduct.formated <- ifelse(any(x.exist),
                              MSCC::adduct.table$Adduct[x.exist],NA)

    data.frame(warning = sum(x.exist )!= 1,
         adduct.input = x,
         adduct.formated = adduct.formated )

  }



  adduct.temp <- sapply(unique(adduct.to.check) ,
              .check_adduct,simplify = F)%>%
    do.call(what = "rbind")

  to.return <-adduct.temp[match(adduct.to.check,adduct.temp$adduct.input),]%>%
    `rownames<-`(NULL)

return(to.return)

}




