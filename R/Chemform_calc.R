#' @describeIn Chemform_calculation chemform_sum
#' @title Chemform Calculation
#' @param ... chemform list or vector
#'
#' @returns chemform
#' @export
#'
#' @examples chemform_sum(chem_formula_template)
chemform_sum <- function(...,return =  c("chemform","matrix")){
  return <- match.arg(return)
  chemform.list <- list(...)
  chemform.list <- unlist(chemform.list,recursive = T)
  chemfomr.m <- chemform_parse(chemform.list,return = "matrix")
  chemfomr.m <- apply(chemfomr.m,2,sum)
  if (return=="chemform") {
   return( chemform_from_ele_matrix(chemfomr.m))
  }

  return(chemfomr.m)
}


#' @describeIn Chemform_calculation chemform_multi
#' @title Chemform Calculation
#' @param chemform chemform
#' @param multi number to multi
#'
#' @returns chemform
#' @export
#'
#' @examples chemform_multi("C6H12O6",1)
chemform_multi <- function(chemform = MSCC::chem_formula_template,
                           multi = 1,
                           return = c("matrix","chemform")){
  return <- match.arg(return)
  if (length(multi)==1 ) {
    multi <- rep(multi,length(chemform))
  }
  chemform.matrix <- chemform_parse(chemform)
  chemform_matrix_multi(chemform.matrix , multi = multi,return = return)


}
