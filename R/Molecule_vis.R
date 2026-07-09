#' Convert SDF to igraph Object
#'
#' @description Converts a chemical structure from SDF format to an igraph object for network visualization.
#'
#' @param sdf An SDF or SDFset object.
#' @param addH Logical, whether to include hydrogen atoms (default FALSE).
#'
#' @return An igraph object representing the molecular structure.
#' @export
get_sdf_igraph <- function(sdf,addH = F){

  .f <- function(sdf,addH){

    #sdf <- canonicalize(sdf)
    #cn <- canonicalNumbering_OB(obmol(sdf))[[1]]
    atom.data <- atomblock(sdf)[,1:2]%>%
      `colnames<-`(c("x","y"))%>%
      as.data.frame()%>%
      rownames_to_column("Atom_id" )#%>%
      #dplyr::mutate(Canonical_Numbering = cn )

    atom.data <- cbind(atom.data,bonds(sdf))%>%
      dplyr::group_by(atom)%>%
      dplyr::mutate(id = Atom_id,
                    label = paste0(" ",atom," "),
                    x = x*100,y = y *100,
                    font.size = 30,
                    borderWidth = 20,
                    font.vadjust = 5,
                    font.align = "center",
                    color.border = "#AAAAAA",
                    color.background = "#FFFFFF",
                    borderWidth = 5,
                    shape = "circle",
                    physics = F)%>%
      dplyr::ungroup()
    bond.data <- bondblock(sdf)[,1:3,drop =F]%>%
      `colnames<-`(c("from","to","bond_type"))%>%
      as.data.frame()%>%
      dplyr::mutate(from = atom.data$Atom_id[from],
                    to = atom.data$Atom_id[to],
                    width = 10*bond_type,
                    color = c("#888888","#444444","#000000")[bond_type],
                    smooth = FALSE
                    )

  #  atom.data <- atom.data%>%
   #   dplyr::arrange(desc(Canonical_Numbering))%>%
   #   dplyr::mutate(RCN_id = paste0(atom,num2str(1:n())),
   #                 label = RCN_id)%>%
   #   dplyr::ungroup()

    sdf.igraph <- igraph::graph_from_data_frame(
      bond.data,vertices = atom.data
    )
    if (!addH) {
      sdf.igraph <- igraph::delete.vertices(sdf.igraph,
                                    V(sdf.igraph)$atom=="H")
    }

    return(sdf.igraph)
  }

  if (class(sdf)=="SDF")
    sdf.igraph <-.f(sdf,addH)


  if (class(sdf)=="SDFset") {
    sdf.igraph <- list()
    sdf.valid <- validSDF(sdf)
    for (i in 1:length(sdf)) {

      sdf.igraph[[i]] <- .f(sdf[[i]],addH  )

    }
  }

  return(sdf.igraph)


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

  sdf.igraph <- sdf.igraph%>%
    sdf_igraph_show_id(show_id)

  vda <- vdata(sdf.igraph)%>%
    dplyr::mutate(x = x-mean(x),
                  y = y-mean(y))%>%
    dplyr::select("id",everything())
  eda <- edata(sdf.igraph)%>%
    dplyr::select(c("from","to","bond_type","width","color","smooth"))

  visNetwork(nodes = vda,edges = eda) %>%
    visLayout( hierarchical = FALSE)

}

vis_cfm_data_fragment <- function(cfmd,fragment_id,
                                  show_id= T,...){

  ig <- cfmd@fragment_igraph[[fragment_id]]%>%
    sdf_igraph_show_id(show_id )
  vis_sdf_igraph(ig,show_id = show_id,...)

}

vis_cfm_data_FG_mapped <- function(cfmd,FG,
                                   show_id= T){

  atom.prob <- get_MSIPAtomMap_fragment_group_atom_map(cfmd,FG)
  get_MSIPAtomMap_sdf_igraph(cfmd)%>%
    sdf_igraph_add_border_color(value = atom.prob)%>%
    vis_sdf_igraph(show_id)
}

vis_cfm_data_fragment_atom_map<- function(cfmd,fragment_id,
                                          show_id= T){

  maps <- cfmd@fragment_atom_map[[fragment_id]]

  sdf.igraphA <- get_cfm_data_sdf_igraph(cfmd,1)%>%
    sdf_igraph_add_border_color(value = rowSums(maps))%>%
    sdf_igraph_show_id(show_id)
  sdf.igraphB <- get_cfm_data_sdf_igraph(cfmd,fragment_id)%>%
    sdf_igraph_add_border_color(value = colSums(maps))%>%
    sdf_igraph_show_id(show_id)

  sdf_igraph_merge(sdf.igraphA,sdf.igraphB)%>%
    vis_sdf_igraph(show_id = show_id)

}

vis_cfm_data_trans_map<- function(cfmd,trans_id,
                                          show_id= T){

  maps <-get_CFM_data_trans_map(cfmd,trans_id )

  fragment.trans <- cfmd@fragment_transition
  sdf.igraphA <- get_cfm_data_sdf_igraph(cfmd,fragment.trans$from[trans_id])%>%
    sdf_igraph_add_border_color(value = rowSums(maps))%>%
    sdf_igraph_show_id(show_id)
  sdf.igraphB <- get_cfm_data_sdf_igraph(cfmd,fragment.trans$to[trans_id])%>%
    sdf_igraph_add_border_color(value = colSums(maps))%>%
    sdf_igraph_show_id(show_id)

  sdf_igraph_merge(sdf.igraphA,sdf.igraphB)%>%
    vis_sdf_igraph(show_id = show_id)

}

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


sdf_igraph_show_id <- function(sdf.igraph,show_id){

  #message_with_time(show_id)
  if (show_id) {
    vdata(sdf.igraph)$label <-   vdata(sdf.igraph)$id%>%
      str_format_len()
    vdata(sdf.igraph)$font.size <- vdata(sdf.igraph)$font.size * 0.5

  }else{
    vdata(sdf.igraph)$label <-  paste0(" ",vdata(sdf.igraph)$atom," ")
  }

  return(sdf.igraph)
}

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



vis_cfm_data_trans_net <- function(cfmd){

  cfmd.trans.ig <- get_CFM_data_trans_igraph(cfmd)

  visNetwork(nodes = vdata(cfmd.trans.ig),
             edges = edata(cfmd.trans.ig))%>%
    visEdges(smooth = F)


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
