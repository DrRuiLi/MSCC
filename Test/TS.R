

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


MSDEV










rt.v <- isotope.matched$feature.rt

rt.dist <- dist(rt.v)

rt.clust <- hclust(rt.dist)
plot(rt.clust)
rect.hclust(rt.clust, h = 5)

a <- isotope.matched%>%
  dplyr::mutate(x = cutree(hclust(dist( feature.rt )),h = 10)   )%>%
  dplyr::distinct(feature.id,.keep_all = T)%>%
  dplyr::arrange(x,-abundance)


msdev <- load_as_var("d:/MSCC.test/MSdev_2023_05_04.Rdata")
compound.table <- readxl::read_excel("d:/MSCC.test/MSCC.compound.xlsx")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(Chem_formula = chemform_formate(Chem_formula),
                Chem_formula_adduct = chemform_adduct(Chem_formula , Adduct)
  )
xcms.xcms <- msdev@xcmsData$positiveMS1
isotopes.table <- chemform_isotopes_pattern_enviPat(compound.table$Chem_formula_adduct[3])
isotopes.table <- match_isotopes_to_featuredef(isotopes.table,
                                                 msdev@xcmsData$positiveFeature)
isotopes.calced <- match_isotopes_to_featureval(isotopes.table,
                                                featureValues(msdev@xcmsData$positiveMS1))

"C23H39D7NO7P" -> chemform
chemform <- chemform_formate(chemform)
chemform
gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)



isotopes.table %>%
  dplyr::filter()


a <- isotope.candidate %>%
  MSdev:::add_multi_column(unique(isotopes.matched$rt.cluster))%>%
  tidyr::pivot_longer(as.character(unique(isotopes.matched$rt.cluster)),names_to = "rt.cluster")%>%
  dplyr::select(-value)%>%
  dplyr::mutate(rt.cluster = as.numeric(rt.cluster))%>%
  dplyr::bind_rows(isotopes.matched)%>%
  dplyr::group_by(rt.cluster,formula)%>%
  dplyr::arrange(rt.cluster,feature.mz)%>%
  dplyr::slice_head(n=1)%>%
  dplyr::ungroup()%>%
  dplyr::arrange(rt.cluster,-abundance)

chemform[!is.na(chemform)] <- chemform_calculate_lc8(Formula1 = chemform[!is.na(chemform)],
                       Formula2 = "H1",
                       sign = 1)%>%
                        as.vector()%>%
                        chemform_formate()
chemform[is.na(chemform)]  <- NA



data(add)




download.file(MONA.URL,destfile = "a.zip",
              method = "wget",)





