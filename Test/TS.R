

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



MS_demo <- load


chemform_isotopes_pattern_enviPat()

system.time(
  enviPat::check_chemform(isotopes = MSCC::isotopes,
                          chemforms = rep("CH2COOH",100000))->a
  )


system.time(
  chemform_parse.for(chemform = rep(chem_formula_template,100))->a
)


system.time(
  chemform_parse.sapple(chemform = rep(chem_formula_template,100))->a
)


system.time(
  chemform_parse.for(chemform = rep(chem_formula_template,1000))->a
)


system.time(
  chemform_formate(chemform = rep(chem_formula_template,10000))->a
)




chemform <- gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
  gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)




grep(pattern = "[[:alpha:]](?=[0-9^a-z])",x=chemform,perl = T, value = T)

gregexec(pattern = "[[:alpha:]](?=[0-9^a-z])",text =chemform,perl = T)



gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
  gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)

str_extract_all(string = chemform,pattern = "[[:alpha:]]((?=[0-9^a-z])|(?![A-Z]))" )

lc8::my_break_formula()
chemform <- chem_formula_template

chemform.ele <- str_extract_all(string = chemform,pattern = "[[:alpha:]]+((?=[0-9^a-z])|(?![A-Z]))" )


for (i in 1:length(chemform)) {

  elements <- unique(chemform.ele[[i]])
  str_extract_all(string = chemform[i],
                  pattern = paste0("(?<=",elements,")[0-9]+"))



}

### 1. number missing
chemform <-"CH2O"


### 2. duplicated
chemform <-"CH3COOHBr[13]CD[14]Br"


# Wed Jun 14 14:49:22 2023 ------------------------------
chemform <- gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
  gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)

str_extract_all(string = chemform,pattern = "[[A-Z]](?=[0-9^A-Z])|[A-Z][a-z]|\\[[0-9]+\\][A-Z](?=[0-9^A-Z])|\\[[0-9]+\\][A-Z][a-z]")%>%
  `names<-`(chemform)->chemform.ele


str_extract_all(string = chemform[i],
                pattern = paste0("(?<=",elements,")[0-9]+"))


system.time(
  chemform_add_num(chemform = rep(chem_formula_template,1000000))->a
)

###

x <- list()
for (i in 1:length(chemform)) {

  elements <- unique(chemform.ele[[i]])
  elements.exp <- elements
  elements.exp[grepl("\\[",elements)] <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements.exp)],
                                              pattern = "[",
                                              replacement = "\\[",fixed = T)
  elements.exp[grepl("\\[",elements)]  <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements.exp)],
                                               pattern = "]",
                                               replacement = "\\]",fixed = T)

  str_extract_all(string = chemform[i],
                  pattern = paste0("(?<=",elements.exp,")[0-9]+"))%>%
    `names<-`(elements)%>%
    sapply(function(x){sum(as.numeric(x))})->x[[i]]



}



####

.get.ele.num <- function(i ){
  elements <- unique(chemform.ele[[i]])
  elements.exp <- elements
  elements.exp[grepl("\\[",elements)] <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements.exp)],
                                              pattern = "[",
                                              replacement = "\\[",fixed = T)
  elements.exp[grepl("\\[",elements)]  <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements.exp)],
                                               pattern = "]",
                                               replacement = "\\]",fixed = T)

  str_extract_all(string = chemform[i],
                  pattern = paste0("(?<=",elements.exp,")[0-9]+"))%>%
    `names<-`(elements)%>%
    sapply(function(x){sum(as.numeric(x))})->x[[i]]



}
sapply(1:length(chemform),.get.ele.num)


# Wed Jun 14 15:55:43 2023 ------------------------------

paste0( "((?<!\\]",elements.exp[!grepl("\\[",elements)],")(?<=",elements.exp[!grepl("\\[",elements)],"))[0-9]+")

str_extract_all(string = chemform,
                pattern = paste0( "((?<!\\]",elements.exp[!grepl("\\[",elements)],")(?<=",elements.exp[!grepl("\\[",elements)],"))[0-9]+"))





chemform <- hmdb_compound_df$chemical_formula


envipat.result <- enviPat::check_chemform(isotopes,chemforms = chemform)

envipat.result <- chemform_formate("(C14H21N1O11)n1H2O1")


envipat.result <- dplyr::mutate(envipat.result,aaaa = new_formula == MSCC)



grepl(pattern = "[^A-z0-9]",x = "C13H21[13](N1)O10")


chemform <- hmdb_compound_df$chemical_formula

a <- chemform_formate(chemform,return = "data.frame")



b <- cbind(a , chemform_parse(chemform , return = "data.frame"))




##################### mz
chemform <- chemform_formate(chem_formula_template)
a <- chemform_mz_lc8(chemform)
b <- chemform_mz(chemform)


ccc <- data.frame(a,b)%>%
  mutate(id = a==b,
         c = abs(a-b))


chemform <- chemform_formate(chem_formula_template)
system.time(
  chemform_mz_lc8(chemform = rep(chemform,100))->a
)
system.time(
  chemform_mz(chemform = rep(chemform,100))->a
)



system.time(
  a <- data.frame(
    c1 = chem_formula_template[sample(1:107,10000,replace = T)],
    c2 = chem_formula_template[sample(1:107,10000,replace = T)]
  )%>%
    dplyr::mutate(merged = chemform_calc(c1,c2))

)

system.time(
  a <- data.frame(
    c1 = chem_formula_template[sample(1:107,10000,replace = T)],
    c2 = "CH2"
  )%>%
    dplyr::mutate(merged = chemform_calc(c1,c2))

)



lapply(1:length(chemform),.get.ele.num)%>%
  do.call("rbind",.)


chemform_parse(chemform,"matrix")

sapply(x, names) %>% unlist() %>% unique()


chemform2.matrix <- chemform_parse(chemform2)

matrix(rep(chemform2.matrix,length(chemform1)),ncol = ncol(chemform2.matrix),byrow = T,
       dimnames = list(NULL,colnames(chemform2.matrix)))


chemform_calc(chem_formula_template,c(""))


system.time(
  a <- chemform_calc(chem_formula_template,
                               chem_formula_template,
                               ".-")
  )


chemform_calc(
  chem_formula_template[1],
              chem_formula_template[1:5],".-")->a

chemform<- "[13]CH-1HHHC-2HD-"







chemform_formate()



elements.exp[!grepl("\\[",elements)] <- paste0( "((?<!\\]",elements.exp[!grepl("\\[",elements)],")(?<=",elements.exp[!grepl("\\[",elements)],"))[\\-0-9]+")
elements.exp[grepl("\\[",elements)] <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                            pattern = "[",
                                            replacement = "\\[",fixed = T)
elements.exp[grepl("\\[",elements)]  <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                             pattern = "]",
                                             replacement = "\\]",fixed = T)
elements.exp[grepl("\\[",elements)]  <- paste0("(?<=",elements.exp[grepl("\\[",elements)] ,")[\\-0-9]+")
# elements.exp[!repl("\\[",elements)] <- paste0("(?<!)")

str_extract_all(string = chemform[i],
                pattern =elements.exp )%>%
  `names<-`(elements)->ele.num




a <- data.frame(cf = chem_formula_template,
                  adduct = sample(adduct.table$Adduct,replace = T,size = length(chem_formula_template)))%>%
  dplyr::mutate(chemform_adduct(chemform = cf,
                                adduct = adduct))


