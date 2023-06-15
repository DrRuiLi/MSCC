

chemform_add_num <- function(chemform){

  ### add number 1 after a char follow:
  #### 1. An alpha, not any lowercase or number exist after it
  ### replace D with [H]

  chemform <- gsub(pattern = "([[:alpha:]](?![0-9^a-z]))" ,replacement = "\\11",x=chemform,perl = T)%>%
    gsub(pattern = "D(?=[0-9])",replacement = "[2]H",perl = T)
  return(chemform)


}

#' chemform_get_ele
#'
#' char considered as a element should follow one of:
#' 1. uppercase with uppercase or number after it
#' 2. continuous uppercase and lowercase
#' 3. continuous \[n] and uppercase with uppercase or number after it
#' 4. continuous \[n] and uppercase and lowercase
#'
#' @param chemform
#'
#' @return
#' @export
#'
#' @examples
#' chemform_get_ele("CCC[13]C")
chemform_get_ele <- function(chemform){



  chemform <- chemform_add_num(chemform)


  chemform.ele <- str_extract_all(string = chemform,pattern = "[[A-Z]](?=[0-9^A-Z])|[A-Z][a-z]|\\[[0-9]+\\][A-Z](?=[0-9^A-Z])|\\[[0-9]+\\][A-Z][a-z]")%>%
    `names<-`(chemform)
  return(chemform.ele)
}



#' chemform_parse
#'
#' parse given chemform to a vector of elements and count
#'
#'
#' @param chemform
#' @param return
#'
#' @return
#' @export
#'
#' @examples
chemform_parse <- function(chemform = chem_formula_template,return = "list"){


  chemform <- chemform_add_num(chemform)
  chemform.ele <- chemform_get_ele(chemform)



  .get.ele.num <- function(i ){

    elements <- unique(chemform.ele[[i]])
    elements.exp <- elements
    elements.exp[!grepl("\\[",elements)] <- paste0( "((?<!\\]",elements.exp[!grepl("\\[",elements)],")(?<=",elements.exp[!grepl("\\[",elements)],"))[0-9]+")
    elements.exp[grepl("\\[",elements)] <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                                pattern = "[",
                                                replacement = "\\[",fixed = T)
    elements.exp[grepl("\\[",elements)]  <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                                 pattern = "]",
                                                 replacement = "\\]",fixed = T)
    elements.exp[grepl("\\[",elements)]  <- paste0("(?<=",elements.exp[grepl("\\[",elements)] ,")[0-9]+")
   # elements.exp[!repl("\\[",elements)] <- paste0("(?<!)")

    str_extract_all(string = chemform[i],
                    pattern =elements.exp )%>%
      `names<-`(elements)->ele.num

    for (j in 1:length(ele.num)) {
      ele.num[[j]] <- sum(as.numeric(ele.num[[j]]))
    }
    unlist(ele.num)


  }



  if (length(chemform) >1) {
    chemform.ele <- sapply(1:length(chemform),.get.ele.num)%>%
      `names<-`(chemform)

  }else{
    chemform.ele <- .get.ele.num(1)
    chemform.ele.matrix <-matrix(chemform.ele,nrow = 1,dimnames = list(chemform,names(chemform.ele)))
    return(chemform.ele.matrix)
  }

  if (return== "matrix") {
    chemform.ele.matrix <- MSdev::list2df(chemform.ele)%>%as.matrix()

    class(chemform.ele.matrix) <- "numeric"
    return(chemform.ele.matrix)
  }

}


#' @title chemform_formate
#' @description
#' formate a chemical formula
#'
#' @param chemform
#'
#' @return
#' @export
#'
#' @examples
chemform_formate <- function(chemform = chem_formula_template,
                             return = "chemform"){

  chemform.raw <- chemform

  chemform <- chemform_add_num(chemform)

  ### is char valid
  char.valid <- !grepl(pattern =  "[^A-z0-9]",x = chemform)


  ### is element valid
  chemform.ele <- chemform_get_ele(chemform)
  ele.all <- MSCC::elem_table$element
  ele.valid <- sapply(chemform.ele,function(x){
    if (is_empty(x)) {
      return(FALSE)
    }
    all(x%in% ele.all)
  })


  ### remove non valid
  chemform[!ele.valid | !char.valid] <-NA


  if (return== "chemform") {

    return(chemform)
  }else{


    return(data.frame(chemform=chemform.raw,
                      formated = chemform,
                      ele.valid ,char.valid))

  }

}




chemform_mz <- function(chemform = chem_formula_template){

  chemform <- chemform_formate(chemform)
  chemform.matrix <- chemform_parse(chemform,return = "matrix")

  ele.mass <- MSCC::elem_table$mass[match(colnames(chemform.matrix),MSCC::elem_table$element)]%>%
    `names<-`(colnames(chemform.matrix))

  mass.matrix <- t(t(chemform.matrix) * ele.mass)
  mass.matrix[is.na(mass.matrix)] <- 0
  chemform.mz <- apply(mass.matrix,1 ,sum )
  return(chemform.mz)

}
