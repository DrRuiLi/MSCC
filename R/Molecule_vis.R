#' Molecule graph from a SMILES string
#'
#' @description Builds an igraph from the canonical RDKit molecule.
#'   Hydrogen atoms are not included.
#'
#' @param sdf A SMILES string, or several SMILES strings.
#' @param addH Ignored. Hydrogens are not part of the graph.
#'
#' @return An igraph object, or a list of them when `sdf` has several strings.
#' @export
get_sdf_igraph <- function(sdf, addH = FALSE) {
  if (!is.character(sdf)) {
    stop("Pass a SMILES string to get_sdf_igraph().", call. = FALSE)
  }
  .one <- function(smiles) {
    g <- .rdkit_parts_igraph(rdkit_molecule_parts(smiles))
    vdata(g) <- vdata(g) %>%
      dplyr::mutate(
        x = x * 100,
        y = y * 100,
        font.size = 30,
        borderWidth = 5,
        font.vadjust = 5,
        font.align = "center",
        color.border = "#AAAAAA",
        color.background = "#FFFFFF",
        shape = "circle",
        physics = FALSE
      )
    ed <- edata(g)
    if (nrow(ed)) {
      bt <- suppressWarnings(as.integer(round(ed$bond_type)))
      bt[!is.finite(bt) | bt < 1L | bt > 3L] <- 1L
      ed$width <- 10 * ed$bond_type
      ed$color <- c("#888888", "#444444", "#000000")[bt]
      ed$smooth <- FALSE
      edata(g) <- ed
    }
    g
  }
  if (length(sdf) == 1L) {
    return(.one(sdf))
  }
  lapply(sdf, .one)
}



#' Visualize SDF igraph
#'
#' @description Creates an interactive visualization of a molecular igraph using visNetwork.
#'
#' @param sdf.igraph An igraph object representing the molecule.
#' @param show_id Logical, whether to show atom IDs as labels (default FALSE).
#' @param ... Additional arguments passed to visNetwork.
#'
#' @return A visNetwork html widget.
#' @export
vis_sdf_igraph <- function(sdf.igraph,show_id = F,...){

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("Package visNetwork is required for vis_sdf_igraph().", call. = FALSE)
  }

  sdf.igraph <- sdf.igraph%>%
    sdf_igraph_show_id(show_id)

  vda <- vdata(sdf.igraph)%>%
    dplyr::mutate(
      x = (x - mean(x, na.rm = TRUE)) * 100,
      y = (y - mean(y, na.rm = TRUE)) * 100,
      physics = FALSE,
      shape = "circle",
      borderWidth = 4,
      label = trimws(as.character(label)),
      font.size = 20,
      font.align = "center",
      font.vadjust = 4,
      font.color = "#1b1b1b"
    )%>%
    dplyr::select("id", dplyr::everything())
  eda <- edata(sdf.igraph)
  if (!"bond_type" %in% names(eda)) eda$bond_type <- 1
  if (!"width" %in% names(eda)) {
    bt <- suppressWarnings(as.numeric(eda$bond_type))
    bt[!is.finite(bt)] <- 1
    eda$width <- 4 * bt
  }
  if (!"smooth" %in% names(eda)) eda$smooth <- FALSE
  if (!"color" %in% names(eda)) eda$color <- "#848484"
  eda <- dplyr::select(
    eda,
    dplyr::any_of(c("from", "to", "bond_type", "width", "color", "smooth"))
  )

  visNetwork::visNetwork(nodes = vda, edges = eda) %>%
    visNetwork::visNodes(
      shape = "circle",
      borderWidth = 4,
      # Equal top/bottom margins grow the circle around the letter.
      # A slightly larger left margin cancels the viewer's leftward bias.
      margin = list(top = 16, right = 12, bottom = 16, left = 18),
      font = list(
        size = 20,
        face = "arial",
        align = "center",
        vadjust = 4,
        color = "#1b1b1b"
      )
    ) %>%
    visNetwork::visPhysics(enabled = FALSE) %>%
    visNetwork::visLayout(hierarchical = FALSE)

}

#' @export
sdf_igraph_add_border_color <- function(sdf.igraph,
                                        value,
                                        color.ramp = colramp(breaks = c(0,Inf,1),
                                                             colors = c("#CCCCCC","#97C2FC","#2B7CE9"))
                                        ){

  sdf.igraph.temp <-sdf.igraph
  ele <-  MSCC::get_sdf_igraph_atom(sdf.igraph)
  col.border <- .get_vis_col(sdf.igraph,value,color.ramp)
  vdata(sdf.igraph.temp)$color.border <- col.border
  return(sdf.igraph.temp)


}
#' Color atom backgrounds on a structure graph.
#'
#' @param sdf.igraph An igraph from a molecular structure.
#' @param value Named numeric weights. Names are atom ids.
#' @param color.ramp Color ramp used by the structure viewer.
#' @return The graph with `color.background` set on vertices.
#' @export
sdf_igraph_add_background_color <- function(sdf.igraph,
                                        value,
                                        color.ramp = colramp(breaks = c(0,Inf,1),
                                                             colors = c("#FFFFFF","#F7844F","#B20C26"))
){

  sdf.igraph.temp <-sdf.igraph
  ele <-  MSCC::get_sdf_igraph_atom(sdf.igraph)
  col.background<- .get_vis_col(sdf.igraph,value,color.ramp)
  vdata(sdf.igraph.temp)$color.background <- col.background
  return(sdf.igraph.temp)


}


#' @export
sdf_igraph_show_id <- function(sdf.igraph, show_id) {
  vd <- vdata(sdf.igraph)
  if (!nrow(vd)) return(sdf.igraph)
  if (isTRUE(show_id)) {
    ids <- vd$id
    if (is.null(ids) || length(ids) != nrow(vd)) ids <- vd$name
    ids <- as.character(ids)
    if (length(ids) != nrow(vd)) ids <- as.character(seq_len(nrow(vd)))
    vd$label <- ids
    # Graphs built from SMILES have no font.size. Multiplying NULL is an
    # empty vector and `$<-.data.frame` then reports "replacement has 0 rows".
    if (!is.null(vd$font.size) && length(vd$font.size) == nrow(vd)) {
      vd$font.size <- vd$font.size * 0.5
    }
  } else if ("atom" %in% names(vd)) {
    vd$label <- trimws(as.character(vd$atom))
  }
  vdata(sdf.igraph) <- vd
  sdf.igraph
}

#' @export
sdf_igraph_merge <- function(sdf.igraphA,sdf.igraphB){


  xa <-vdata(sdf.igraphA)$x
  xb <- vdata(sdf.igraphB)$x

  vdata(sdf.igraphA)$x <- xa-max(xa)-diff(range(xa))/5
  vdata(sdf.igraphB)$x <- xb- min(xb)+diff(range(xb))/5

  vdata(sdf.igraphA)$name <- paste0("F_",vdata(sdf.igraphA)$name)
  vdata(sdf.igraphB)$name <- paste0("T_",vdata(sdf.igraphB)$name)

  vdata(sdf.igraphA)$id <- paste0("F_",vdata(sdf.igraphA)$id)
  vdata(sdf.igraphB)$id <- paste0("T_",vdata(sdf.igraphB)$id)

  edata(sdf.igraphA)$id <- paste0("F_",edata(sdf.igraphA)$id)
  edata(sdf.igraphB)$id <- paste0("T_",edata(sdf.igraphB)$id)

  nodes <- bind_rows(
    vdata(sdf.igraphA),
    vdata(sdf.igraphB)
  )
  #nodes$size <- 20
  edges <- rbind(
    edata(sdf.igraphA),
    edata(sdf.igraphB)
  )
  new.igraph <- igraph::graph_from_data_frame(edges,
                                      vertices = nodes)
  return(new.igraph)
}


#' Place several structure graphs side by side.
#'
#' @param ... Structure graphs.
#' @param vertical Stack along y when `TRUE`, otherwise along x.
#' @return One igraph containing every input structure.
#' @export
sdf_igraph_merge_all <- function(...,vertical = T){


  ### id for input mig
  {
    ig.list <- c(...)
    if (is.null(names(ig.list))) {
      names(ig.list) <- paste0("CP",1:length(ig.list))
    }
  }

  ### new id and name
  {
    ig.list <- sapply(names(ig.list),function(id){
      ig <- ig.list[[id]]
      vdata(ig)$name <- paste0(id,"_",vdata(ig)$name)
      vdata(ig)$id <- paste0(id,"_",vdata(ig)$id)
      edata(ig)$id <- paste0(id,"_",edata(ig)$id)
      ig
    }, USE.NAMES =T)
  }


  ### new position
  {
    ig.size <- sapply(names(ig.list),function(id){

      ig <- ig.list[[id]]
      s <- c(x = diff(range(vdata(ig)$x)),
                 y = diff(range(vdata(ig)$y) ))
      if (vertical) {
        return(unname(s["y"]))
      }else{
        return(unname(s["x"]))
      }})

    for (i in seq_along(ig.list)) {

      if (vertical) {

        vdata(ig.list[[i]])$x <- vdata(ig.list[[i]])$x - mean(vdata(ig.list[[i]])$x)
        vdata(ig.list[[i]])$y <- vdata(ig.list[[i]])$y - mean(vdata(ig.list[[i]])$y)
        vdata(ig.list[[i]])$y <- vdata(ig.list[[i]])$y +cumsum(ig.size)[i]

      }else{

        vdata(ig.list[[i]])$x <- vdata(ig.list[[i]])$x - mean(vdata(ig.list[[i]])$x)
        vdata(ig.list[[i]])$y <- vdata(ig.list[[i]])$y - mean(vdata(ig.list[[i]])$y)
        vdata(ig.list[[i]])$x <- vdata(ig.list[[i]])$x +cumsum(ig.size)[i]
      }
    }



  }


  ### combine
  {
    nodes <- lapply(ig.list,vdata)%>%
      do.call(bind_rows,.)
    edges <- lapply(ig.list,edata)%>%
      do.call(bind_rows,.)

    new.igraph <- igraph::graph_from_data_frame(edges,
                                                vertices = nodes)
  }

  return(new.igraph)


}



.get_vis_col <- function(sdf.igraph,prob,colramp = colramp(),na.col = "#AAAAAA"){

  ele <- MSCC::get_sdf_igraph_atom(sdf.igraph,ele = "all")
  if (is.null(prob)) prob <- 0
  if (is.null(names(prob))&length(prob)==1){
    prob <- rep(prob,length(ele))
    names(prob)<-ele
  }else{
    x <- setdiff(ele,names(prob))
    xx <- rep(0,length(x))
    names(xx) <-x
    prob <- c(prob,xx)[ele]
  }
  col <- colramp(prob)
  names(col)<-names(prob)
  col[is.na(col)] <- na.col
  col
}



rotate_nodes_safe <- function(ig, angle_deg) {

  theta <- angle_deg * (pi / 180)

  nodes <- vdata(ig)
  range_x <- max(nodes$x) - min(nodes$x)
  range_y <- max(nodes$y) - min(nodes$y)
  avg_range <- (range_x + range_y) / 2

  cx <- mean(nodes$x)
  cy <- mean(nodes$y)
  nodes$x_std <- (nodes$x - cx)
  nodes$y_std <- (nodes$y - cy)

  nodes$x_new <- nodes$x_std * cos(theta) - nodes$y_std * sin(theta)
  nodes$y_new <- nodes$x_std * sin(theta) + nodes$y_std * cos(theta)

  nodes$x <- nodes$x_new + cx
  nodes$y <- nodes$y_new + cy

  nodes -> vdata(ig)
  return(ig)
}
