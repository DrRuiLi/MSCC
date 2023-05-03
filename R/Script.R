

ff <- function(xx){

  aa <- elem_table[,1:5]%>%
    as.data.frame()

  bb <- MSdev::edit_df_in_excel(isotopes)

  enviPat::isopattern(isotopes = aa ,chemforms = "C22H40O12",

                      threshold=0.1,
                      plotit=TRUE,
                      charge=FALSE,
                      emass=0.00054858,
                      algo=1)
  enviPat::isopattern(isotopes ,chemforms = "C22H40O12D1",

                      threshold=0.1,
                      plotit=TRUE,
                      charge=FALSE,
                      emass=0.00054858,
                      algo=1)

  iso.table <-elem_table[,1:5]%>%as.data.frame()
  use_data(iso.table)

  enviPat::isopattern(isotopes = isotopes,
                      chemforms = chemform,
                      threshold = 0.01 )[[1]]




  compound.template <- hmdb_compound_df%>%
    dplyr::filter(!is.na(kegg_id))%>%
    dplyr::select(Compound = name,
                  Chem_formula = chemical_formula)%>%
    dplyr::mutate(.after = Chem_formula,
                  Adduct = "M+H",
                  Iso_label = NA)

  a <- check_chemform(isotopes ,compound.template$Chem_formula)

  isotopes_pattern_enviPat(a)

  a <- readxl::read_excel("d:/isoadj.test/isoadj.data.xlsx")


  iso.data <- readxl::read_excel("d:/isoadj.test/isoadj.data.xlsx")
  iso.compound <- readxl::read_excel("d:/isoadj.test/isoadj.compound.xlsx")



}
