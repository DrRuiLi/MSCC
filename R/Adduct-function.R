### get adduct.table from package enviPat, and add/correct some info
### adduct.table will be stored in /data
### a column "Adduct_Syn" for char match, only string list in Adduct_Syn will be considered as valid adduct
get.adduct.table.from.enviPat <- function(){


  data("adducts",package = "enviPat")

  adduct.table <- adducts%>%
    dplyr::add_row(Name = "M-H2O+H",calc = "M-14.987089588",
                   Charge=1,Mult=1,Mass=chemform_mz_lc8("H-1O-1"),Ion_mode="positive",
                   Formula_add="FALSE",Formula_ded="H1O1",Multi=1)%>%
    dplyr::add_row(Name = "M+H2O+H",calc = "M+19.01838972",
                   Charge=1,Mult=1,Mass=chemform_mz_lc8("H3O1"),Ion_mode="positive",
                   Formula_add="H3O1",Formula_ded="FALSE",Multi=1)%>%
    dplyr::filter(Name != "2M+3H2O+2H")%>%
    dplyr::rowwise()%>%
    dplyr::mutate(a=case_when(Formula_add=="FALSE"~"C0",
                              T~ Formula_add),
                  b=case_when(Formula_ded=="FALSE"~"C0",
                              T~ Formula_ded),
                  c = case_when(abs(Charge) == 1~ "",
                                T~ as.character(abs(Charge))),
                  Formula_diff = chemform_calc(
                    chemform_formate(a),
                    chemform_formate(b),
                    "-",return = "chemform"
                  ),
                  Adduct = paste0("[",Name,"]",c,ifelse(Ion_mode=="positive","+","-"))
    )%>%
    dplyr::select(-a,-b,-c,-Mult)

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
                    ))%>%
      dplyr::mutate(Formula_diff = case_when(Formula_diff==""~"C0",
                                             T~Formula_diff))

    }

  #  use_data(adduct.table,overwrite = T)

  return(adduct.table)
}






#' @title chemform_adduct_check
#' @description
#' check if string is a adduct (match in adduct.table)
#'
#'
#' @param adduct.to.check
#'
#' @return df
#' @export
#' @import tidyverse

chemform_adduct_check <- function(adduct.to.check ){

  .check_adduct <- function(x){

    grepl(x = MSCC::adduct.table$Adduct_Syn,
          pattern = paste0(x,";"),
          fixed = T)->x.exist
    adduct.formated <- ifelse(any(x.exist),
                              MSCC::adduct.table$Adduct[x.exist],NA)

    data.frame(warning = sum(x.exist )!= 1,
               adduct.input = x,
               adduct.formated = adduct.formated )%>%
      dplyr::mutate(MSCC::adduct.table[match(adduct.formated,MSCC::adduct.table$Adduct),]  )

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
#' @import tidyverse
#' @return mz or df
#' @export
#'

chemform_adduct <- function(chemform = chem_formula_template,
                            adduct = "[M+H]+",
                            value = c("mz","chemfrom","all")){

  value <- match.arg(value)
  ### formate and unique
  {
    chemfrom_raw <- chemform
    chemform <- chemform_formate(chemform)
    chemform.f <- factor(chemform)

  }


  ###construct df
  {

    adduct.check <- chemform_adduct_check(adduct)
    adduct.check <- adduct.check[!adduct.check$warning,]

    chem_df <- expand.grid(id = 1:length(chemform),
                           adduct = 1:nrow(adduct.check),
                           stringsAsFactors = F)%>%
      dplyr::mutate(chemform = chemform[id],
                    adduct.check[adduct,c("Formula_diff","Multi","Charge")]
                    )

  }

  ### matrix calculation
  {
    chemform.matrix <- chemform_parse(chemform)
    adduct.matrix <- chemform_parse(adduct.check$Formula_diff)

    chemform.matrix <- chemform.matrix[chem_df$id ,,drop=F ]
    adduct.matrix <- adduct.matrix[chem_df$adduct,,drop=F ]
    chemform.matrix.multi <- chemform_matrix_multi(chemform.matrix,
                                                   chem_df$Multi,return = "matrix")
    chemform.matrix.calc <- chemform_matrix_calc(chemform.matrix.multi,
                                                 adduct.matrix,
                                                 calc = "+",
                                                 return = "matrix")

    chemform.matrix.mz <- chemform_matrix_mz(chemform.matrix.calc,
                                             charge = chem_df$Charge)
  }

  ### remove and return
  {

    chem_df$chemform.adduct <- chemform_from_ele_matrix(chemform.matrix.calc)
    chem_df$chemform.adduct.mz <- chemform.matrix.mz
    chem_df$chemform.raw <- chemfrom_raw[chem_df$id]
    chem_df$adduct <- adduct.check$Adduct[chem_df$adduct]
    ele.count.valid <- apply(chemform.matrix.calc, 1, function(x){
      all(x>=0)
    } )
    idx.error <- !ele.count.valid|is.na(chem_df$chemform)
    chem_df$chemform.adduct.mz[idx.error] <- NA
    chem_df$chemform.adduct[idx.error] <- NA

  }


  ### rerturn
  {
    to.return <- switch (value,
                         "all" = chem_df[!idx.error,],
                         "mz" = chem_df$chemform.adduct.mz,
                         "chemform" = chem_df$chemform.adduct
    )

    return(to.return)



  }




}






