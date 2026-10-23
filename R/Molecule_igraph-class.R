setClass("Molecule_igraph",
         slots = list(
           molecule_info = "list",
           igraph = "ANY",
           isotopomer = "data.frame"
         ))

#' @export
Molecule_igraph <- function(){
  new("Molecule_igraph")
}


setMethod("show",
  "Molecule_igraph",
  definition = function(object) {
    mf <- object@molecule_info$formula
    if (is.null(mf) || !length(mf) || is.na(mf[[1]]) || !nzchar(mf[[1]])) {
      smi <- object@molecule_info$smiles
      mf <- if (is.null(smi) || !length(smi) || is.na(smi[[1]])) "unknown" else smi[[1]]
    }
    cat(
      "Molecule_igraph: ", paste(mf, collapse = ", "), " ",
      nrow(object@isotopomer), " isotopomers\n",
      sep = ""
    )
    invisible(object)
  }
)




setMethod("vdata",
  "Molecule_igraph",
  definition = function(object) {
    vdata(object@igraph)
  }
)




setMethod("vdata<-",
  "Molecule_igraph",
  definition = function(object, value) {
    vdata(object@igraph) <- value
    object
  }
)



setMethod("edata",
  "Molecule_igraph",
  definition = function(object) {
    edata(object@igraph)
  }
)



setMethod("edata<-",
  "Molecule_igraph",
  definition = function(object, value) {
    edata(object@igraph) <- value
    object
  }
)





setMethod(f = "plot",
  signature = "Molecule_igraph",
  definition = function(object,x,y) {
    smi <- object@molecule_info$smiles
    if (is.null(smi) || !length(smi) || is.na(smi[[1]]) || !nzchar(smi[[1]])) {
      stop("This molecule has no SMILES to plot.", call. = FALSE)
    }
    print(rdkit_plot_mol(rdkit_mol_from_smiles(smi[[1]])))
    invisible()
  }
)



setMethod("atom",
          "Molecule_igraph",
          definition = function(object,
                                element = "ANY"){

            if (element=="ANY") {
              data(element_table)
              element <- element_table$element
            }
            vdata(object)%>%
              dplyr::filter(element %in% !!element)%>%
              dplyr::pull(name)
          })


setMethod("get_element",
          "Molecule_igraph",
          definition = function(object,...){
            vdata(object)%>%
              dplyr::pull(element)
          })


setMethod("formula","Molecule_igraph",
          definition = function(x,...){
            mf <- x@molecule_info$formula
            if (is.null(mf) || !length(mf)) NA_character_ else as.character(mf[[1]])
          })


setMethod("as.character",
          "Molecule_igraph",
          definition =  function(x,...){
            print(paste0("Molecule_igraph: ", formula(x)))
          })



