#' Vertex attribute table for graph-like objects
#'
#' @param object Graph-like object (for example an `igraph`).
#'
#' @return Data frame of vertex attributes.
#' @export
setGeneric("vdata", function(object) igraph::as_data_frame(object, "vertices"))

#' Edge attribute table for graph-like objects
#'
#' @param object Graph-like object (for example an `igraph`).
#'
#' @return Data frame of edge attributes.
#' @export
setGeneric("edata", function(object) igraph::as_data_frame(object, "edges"))

#' Replace vertex attributes for graph-like objects
#'
#' @param object Graph-like object.
#' @param value Data frame-like vertex attributes.
#'
#' @return Updated object.
#' @export
setGeneric("vdata<-", function(object, value) {
  if (inherits(object, "igraph")) {
    igraph::vertex.attributes(object) <- as.list(value)
    return(object)
  }
  standardGeneric("vdata<-")
})

#' Replace edge attributes for graph-like objects
#'
#' @param object Graph-like object.
#' @param value Data frame-like edge attributes.
#'
#' @return Updated object.
#' @export
setGeneric("edata<-", function(object, value) {
  if (inherits(object, "igraph")) {
    value <- value[, !grepl("^from$|^to$", colnames(value)), drop = FALSE]
    igraph::edge.attributes(object) <- as.list(value)
    return(object)
  }
  standardGeneric("edata<-")
})

#' Atom accessor generic
#'
#' @param object Graph-like object.
#' @param element Optional element filter.
#'
#' @return Atom identifiers.
#' @export
setGeneric("atom", function(object, element = "ANY") standardGeneric("atom"))

#' Element accessor generic
#'
#' @param object Graph-like object.
#' @param ... Additional arguments.
#'
#' @return Element vector.
#' @export
setGeneric("get_element", function(object, ...) standardGeneric("get_element"))

#' Create RXNMapper callable
#'
#' @return Python callable from `rxnmapper`.
#' @export
get_RXNMapper <- function() {
  rxnmp <- reticulate::py_suppress_warnings(reticulate::import("rxnmapper"))
  rxn_mapper <- rxnmp$RXNMapper()
  rxn_mapper$get_attention_guided_atom_maps
}
