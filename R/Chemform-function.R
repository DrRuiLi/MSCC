
chemform_add_num <- function(chemform){

  ### add number 1 after a char follow:
  #### 1. An alpha, not any lowercase or number exist after it
  ### replace D with [H]

  chemform <- gsub(pattern = "([[:alpha:]](?![0-9^a-z\\-]))",x= chemform ,replacement = "\\11",perl = T)
  chemform<- gsub(pattern = "[\\-](?![0-9])",x= chemform ,replacement = "\\-1",perl = T)
  chemform<-   gsub(pattern = "D(?=[0-9\\-])",x= chemform ,replacement = "[2]H",perl = T)


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
#' @return list of elements
#' @import magrittr
#' @export
#'


chemform_get_ele <- function(chemform){



  chemform <- chemform_add_num(chemform)


  chemform.ele <- stringr::str_extract_all(string = chemform,pattern = "[[A-Z]](?=[0-9\\-^A-Z])|[A-Z][a-z]|\\[[0-9]+\\][A-Z](?=[0-9\\-^A-Z])|\\[[0-9]+\\][A-Z][a-z]")%>%
    `names<-`(chemform)
  return(chemform.ele)
}



#' chemform_parse
#'
#' parse given chemform to a vector of elements and count
#'
#'
#' @param chemform formula
#' @param return return type
#' @import tidyverse magrittr
#' @return ele matrix of list
#' @export
#'

chemform_parse <- function(chemform = MSCC::chem_formula_template,return = "matrix"){



  chemform <- chemform_formate(chemform)
  ### unique
  if (any(duplicated(chemform))) {
    chemform_uni <- unique(chemform)
    chemform.parsed <- chemform_parse(chemform_uni,return=return)
    chemform[is.na(chemform)] <- "NA"
    if (return== "matrix"){
      return(chemform.parsed[chemform,,drop =F])
  }else{
    return(chemform.parsed[chemform])
    }
  }



  chemform.ele <- chemform_get_ele(chemform)



  .get.ele.num <- function(i ){

    elements <- unique(chemform.ele[[i]])
    if (length(elements)==0|all(is.na(elements))) {
      return(c("C"=0))
    }
    elements.exp <- elements
    elements.exp[!grepl("\\[",elements)] <- paste0( "((?<!\\]",elements.exp[!grepl("\\[",elements)],")(?<=",elements.exp[!grepl("\\[",elements)],"))[0-9\\-]+")
    elements.exp[grepl("\\[",elements)] <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                                pattern = "[",
                                                replacement = "\\[",fixed = T)
    elements.exp[grepl("\\[",elements)]  <- gsub(x = elements.exp[grepl(pattern = "\\[",x = elements)],
                                                 pattern = "]",
                                                 replacement = "\\]",fixed = T)
    elements.exp[grepl("\\[",elements)]  <- paste0("(?<=",elements.exp[grepl("\\[",elements)] ,")[0-9\\-]+")
    # elements.exp[!repl("\\[",elements)] <- paste0("(?<!)")

    stringr::str_extract_all(string = chemform[i],
                             pattern =elements.exp )%>%
      `names<-`(elements)->ele.num

    for (j in 1:length(ele.num)) {
      ele.num[[j]] <- sum(as.numeric(ele.num[[j]]))
    }
    unlist(ele.num)


  }




    chemform.ele.count <-
      lapply(1:length(chemform),.get.ele.num) %>%
      `names<-`(chemform)


  if (return== "matrix") {
    chemform.ele.matrix <- MSdev::list2df(chemform.ele.count)%>%as.matrix()

    class(chemform.ele.matrix) <- "numeric"
    chemform.ele.matrix[is.na(chemform.ele.matrix)] <-0
    return(chemform.ele.matrix)
  }

  return(chemform.ele.count)

}


#' @title chemform_formate
#' @description
#' formate a chemical formula
#'
#' @param chemform
#'
#' @return str or df
#' @export
#'
#' @import tidyverse
chemform_formate <- function(chemform = MSCC::chem_formula_template,
                             return = c("chemform","all")){

  return <- match.arg(return)

  ### unique
  if (any(duplicated(chemform))) {
    chemform_uni <- unique(chemform)
    chemform.formated <- chemform_formate(chemform_uni,return=return)
    if (return== "chemform"){
      chemform.formated <- chemform.formated[match(chemform,chemform_uni)]
      return(chemform.formated)
    }else{
      chemform.formated <- chemform.formated[match(chemform,chemform.formated$chemform),]
      return(chemform.formated)
    }
  }




  ### valid
  {
    chemform.raw <- chemform
    chemform <- chemform_add_num(chemform)

    ### is char valid
    char.valid <- !grepl(pattern =  "[^A-z0-9\\-]",x = chemform)


    ### is element valid
    chemform.ele <- chemform_get_ele(chemform)
    ele.all <- MSCC::elem_table$element
    ele.valid <- sapply(chemform.ele,function(x){
      if (rlang::is_empty(x)) {
        return(FALSE)
      }
      all(x%in% ele.all)
    })%>%unname()

  }



  ### remove non valid
  chemform[!ele.valid | !char.valid] <-NA


  if (return== "chemform") {

    return(chemform)
  }else{


    return(data.frame(chemform = chemform.raw,
                      formated = chemform,
                      ele.valid ,char.valid))

  }

}



chemform_from_ele_matrix <- function(x){

  if (is.null(dim(x))) {
    y<-x
    dim(y) <- c(1,length(x))
    colnames(y) <- names(x)
    x <- y
  }
  apply(x ,1, function(x){
    x <- x[x!=0]
    paste0(paste0(names(x),x),collapse = "")
  }) -> chemform
  return(chemform)

}


#' chemform_mz
#'
#' calculate mz of a given chemform
#'
#' @param chemform
#'
#' @return mz
#' @export
#' @import tidyverse
#'

chemform_mz <- function(chemform = MSCC::chem_formula_template,
                        charge = 0){

  if (length(charge)==1 ) {
    charge <- rep(charge,length(chemform))
  }

  chemform.matrix <- chemform_parse(chemform,return = "matrix")

  chemform.mz <- chemform_matrix_mz(chemform.matrix,charge)
  return(chemform.mz)

}


chemform_matrix_mz <- function(chemform.matrix = chemform_parse(MSCC::chem_formula_template),
                               charge = 0){
  if (length(charge)==1 ) {
    charge <- rep(charge,nrow(chemform.matrix))
  }

  ele.mass <- MSCC::elem_table$mass[match(colnames(chemform.matrix),MSCC::elem_table$element)]%>%
    `names<-`(colnames(chemform.matrix))

  mass.matrix <- t(t(chemform.matrix) * ele.mass)
  mass.matrix[is.na(mass.matrix)] <- 0
  chemform.mz <- apply(mass.matrix,1 ,sum )%>%unname()
  e_mass = 0.00054857990943
  e_diff <- charge
  e_charge <- abs(charge)
  e_charge[e_charge==0] <- 1
  chemform.mz = (chemform.mz - e_mass * e_diff)/e_charge
  return(chemform.mz)


}


#' chemform_calc
#'
#' chemform calculator, calc "+" or "-" for vector operation, which return a vetor with length of m.
#' calc ".+" or ".-" for matrix operation, which return a matrix m x n
#'
#'
#'
#' @param chemform1 vector of chemform, length as m
#' @param chemform2 vector of chemform, length as n
#' @param calc + or -
#' @import tidyverse magrittr
#' @return chemform vector or matrix
#' @export
#'

chemform_calc <- function(chemform1 = MSCC::chem_formula_template ,
                          chemform2 = rev(MSCC::chem_formula_template),
                          calc = "+" ,
                          return = c("matrix","chemform")){
  return <- match.arg(return)
  if (length(chemform1)!=length(chemform2)&length(chemform2)!=1&calc%in% c("+","-")) {
    stop(crayon::yellow("chemform_calc input error: when calc as + or - , length of chemform2 should be 1 or same with chemform1"))

  }
  chemform1.matrix <- chemform_parse(chemform1,return = "matrix")
  if (length(chemform2)==1& calc %in% c("+","-")) {
    chemform2.matrix <- chemform_parse(chemform2)
    chemform2.matrix<-matrix(rep(chemform2.matrix,length(chemform1)),ncol = ncol(chemform2.matrix),byrow = T,
           dimnames = list(NULL,colnames(chemform2.matrix)))
  }else{

    chemform2.matrix <- chemform_parse(chemform2,return = "matrix")

  }

  chemform.calced <- chemform_matrix_calc(chemform1.matrix,
                                          chemform2.matrix,
                                          calc= calc,
                                          return=return)

  return(chemform.calced)


}


### foundation of chemform_calc()
chemform_matrix_calc <- function(chemform1.matrix,
                                 chemform2.matrix,
                                 calc = "+",
                                 return = c("matrix","chemform")
){
  return <- match.arg(return)
  ele.all <- union(colnames(chemform1.matrix),
                   colnames(chemform2.matrix))


  chemform1.matrix <- cbind(chemform1.matrix,
                            matrix(data = 0,nrow = nrow(chemform1.matrix),
                                   ncol = length(setdiff(ele.all,colnames(chemform1.matrix))),
                                   dimnames = list(rownames(chemform1.matrix),
                                                   setdiff(ele.all,colnames(chemform1.matrix))))
  )
  chemform2.matrix <- cbind(chemform2.matrix,
                            matrix(data = 0,nrow = nrow(chemform2.matrix),
                                   ncol = length(setdiff(ele.all,colnames(chemform2.matrix))),
                                   dimnames = list(rownames(chemform2.matrix),
                                                   setdiff(ele.all,colnames(chemform2.matrix))))
  )

  chemform1.matrix <- chemform1.matrix[,ele.all,drop = F]
  chemform2.matrix <- chemform2.matrix[,ele.all,drop = F]

  if (calc == "+") {
    chemform.matrix.calced <- chemform1.matrix+chemform2.matrix
    if (return == "matrix") {
      return(chemform.matrix.calced)
    }
    chemform.calced <- chemform_from_ele_matrix(chemform.matrix.calced)
    names(chemform.calced) <-NULL
  }

  if (calc == "-") {
    chemform.matrix.calced <- chemform1.matrix-chemform2.matrix
    if (return == "matrix") {
      return(chemform.matrix.calced)
    }
    chemform.calced <- chemform_from_ele_matrix(chemform.matrix.calced)
    names(chemform.calced) <-NULL
  }


  if (calc == ".+") {
    chemform1.matrix.expand <-chemform1.matrix[rep(1:nrow(chemform1.matrix),each = nrow(chemform2.matrix)),]
    chemform2.matrix.expand<- chemform2.matrix[rep(1:nrow(chemform2.matrix),times = nrow(chemform1.matrix)),]
    chemform.matrix.expand.calced <- chemform1.matrix.expand+chemform2.matrix.expand
    chemform.expand.calced <- chemform_from_ele_matrix(chemform.matrix.expand.calced)
    chemform.calced <- matrix(chemform.expand.calced,
                              nrow = length(chemform1),byrow = T)
  }

  if (calc == ".-") {
    chemform1.matrix.expand <-chemform1.matrix[rep(1:nrow(chemform1.matrix),each = nrow(chemform2.matrix)),]
    chemform2.matrix.expand<- chemform2.matrix[rep(1:nrow(chemform2.matrix),times = nrow(chemform1.matrix)),]
    chemform.matrix.expand.calced <- chemform1.matrix.expand - chemform2.matrix.expand
    chemform.expand.calced <- chemform_from_ele_matrix(chemform.matrix.expand.calced)
    chemform.calced <- matrix(chemform.expand.calced,
                              nrow = nrow(chemform1.matrix),byrow = T)
  }


  return(chemform.calced)

}




chemform_matrix_multi <- function(chemform.matrix = chemform_parse(MSCC::chem_formula_template),
                                  multi = 1,
                                  return = c("chemform","matrix")){
  return <- match.arg(return)
  if (length(multi)==1 ) {
    multi <- rep(multi,nrow(chemform.matrix))
  }
  chemform.matrix <- chemform.matrix*multi
  if (return == "chemform") {
    return(chemform_from_ele_matrix(chemform.matrix))
  }
  return(chemform.matrix)

}


