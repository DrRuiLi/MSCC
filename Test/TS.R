

gsub("([a])", "\\1_\\1_", "abc and ABC")


grep(x = compound.formula$Chem_formula,
      pattern = "[A-Z][A-Z]",value = T)


sub(x = "C26H45NO7S",
    pattern ="[A-Z][A-Z]",
    replacement = "1")


str <-"C26H45NrO7S"

str_locate_all(string = str,
           pattern = "([A-z][^0-9|a-z])|[A-z]$")

str_locate_all(string = str,
               pattern = "[A-z](?![0-9^a-z])")

str_insrt <- function(str,loc = loc){
  str_new <- paste0(strsplit())
}


data("adducts",package = "enviPat")

c <- adducts%>%
  dplyr::rowwise()%>%
  dplyr::mutate(a=case_when(Formula_add=="FALSE"~"C0",
                                       T~ Formula_add),
                b=case_when(Formula_ded=="FALSE"~"C0",
                                       T~ Formula_ded),
                Formula_diff = chemform_calculate_lc8(
                  chemform_formate(a),
                  chemform_formate(b),
                  -1
                )
               )%>%
  dplyr::select(-a,-b)


xcms.features <- xcms.features%>%
  dplyr::mutate(ppm = (mzmax-mzmin)/mzmed/1e-6)




b <- a%>%
  pivot_longer(11:12)%>%
  dplyr::mutate(feature.mz = str_short(feature.mz,8),
                feature.mz =  factor(feature.mz))

ggplot(b)+
  geom_bar(aes(x = feature.mz,y= value,fill = name),stat = "identity",position = "dodge")+
  labs(x=NULL,y = NULL)+
  theme(axis.text.x = element_text(angle = -45))->p


open_ggplot_win(p,5,3)

