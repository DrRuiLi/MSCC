### get adduct.table from package enviPat, and add/correct some info
### adduct.table will be stored in /data
### a column "Adduct_Syn" for char match, only string list in Adduct_Syn will be considered as valid adduct
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
                            adduct = "[M+H]+",
                            value = "all"){



  chemform <- chemform_formate(chemform)

  if (length(adduct)==1 ) {
    adduct <- rep(adduct,length(chemform))
  }
  #adduct <- sample(MSCC::adduct.table$Adduct,replace = T,length(chemform))


  adduct.check <- chemform_adduct_check(adduct)


  chem_df <- data.frame(chemform = chemform,
                        adduct = adduct.check$adduct.formated)%>%
    dplyr::mutate(chemform = chemform,
                  MSCC::adduct.table[match(adduct,MSCC::adduct.table$Adduct),c("Formula_diff","Multi","Charge")],
                  chemform.adduct = case_when(!is.na(chemform)~chemform_multi(chemform,Multi)),
                  chemform.adduct = case_when(!is.na(chemform)~chemform_calc(chemform.adduct,Formula_diff,calc = "+")),
                  chemform.adduct.mz = case_when(!is.na(chemform)~chemform_mz(chemform.adduct,Charge)),
                  chemform.adduct.mass = case_when(!is.na(chemform)~chemform.adduct.mz * abs(Charge))
    )




  to.return <- switch (value,
    "all" = chem_df,
    "mz" = chem_df$chemform.adduct.mz,
    "chemform" = chem_df$chemform.adduct
  )

  return(to.return)


}


chemform_multi <- function(chemform = chem_formula_template,
                           multi = 1){
  if (length(multi)==1 ) {
    multi <- rep(multi,length(chemform))
  }

  .cm <- function(c,m){
    cf <- ""
    if (m<1|is.na(m)) {
      return(NA)
    }
    for (i in 1:m) {
      cf <- chemform_calc(cf,c,calc = "+")
    }
    return(cf)
  }

  sapply(1:length(chemform),function(i){
    .cm(chemform[i],multi[i])
  })

}




