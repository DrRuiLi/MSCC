# Legacy migration anchor from MSdev `R/0_class.R`.
# Core graph generics/classes are implemented in `graph_primitives.R`.

if (!methods::isGeneric("vdata")) {
  methods::setGeneric("vdata", function(object) standardGeneric("vdata"))
}

if (!methods::isGeneric("edata")) {
  methods::setGeneric("edata", function(object) standardGeneric("edata"))
}

if (!methods::isGeneric("vdata<-")) {
  methods::setGeneric("vdata<-", function(object, value) standardGeneric("vdata<-"))
}

if (!methods::isGeneric("edata<-")) {
  methods::setGeneric("edata<-", function(object, value) standardGeneric("edata<-"))
}

if (!methods::isGeneric("atom")) {
  methods::setGeneric("atom", function(object, element = "ANY") standardGeneric("atom"))
}

if (!methods::isGeneric("get_element")) {
  methods::setGeneric("get_element", function(object, ...) standardGeneric("get_element"))
}
