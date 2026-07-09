setClass("Molecule_igraph",
         slots = list(
           molecule_info = "list",
           sdf = "SDF",
           igraph = "ANY",
           isotopomer = "data.frame"
         ))

#' @importClassesFrom ChemmineR SDF
#' @export
Molecule_igraph <- function(){
  new("Molecule_igraph")
}


setMethod("show",
  "Molecule_igraph",
  definition = function(object) {
    print(paste0("Molecule_igraph: ", unname(MSCC::MF2(object@sdf, addH =  T))," ",
                 nrow(object@isotopomer)," isotopomers"
                 ))

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
    plot(object@sdf)
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
            unname(MF(x@sdf,addH = T))
          })


setMethod("as.character",
          "Molecule_igraph",
          definition =  function(x,...){
            print(paste0("Molecule_igraph: ", unname(MF(x@sdf, addH =  T))))
          })



