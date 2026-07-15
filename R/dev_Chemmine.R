ggplot_sdf <- function(sdf,
                       cex = 1,
                       show_ele = F){


  sdf.formula <- MF(sdf,addH=T)
  sdf.mz <- MSCC::chemform_mz(sdf.formula)%>%round(digits = 4)
  atom.data <- atomblock(sdf)[,1:2]%>%
    `colnames<-`(c("x","y"))%>%
    as.data.frame()%>%
    rownames_to_column("Atom_id" )%>%
    dplyr::mutate(element = str_extract(Atom_id,
                                        "[:alpha:]*"))
  bond.length.short <- ifelse(show_ele,0.1,0)
  bond.data <- bondblock(sdf)[,1:3]%>%
    `colnames<-`(c("from","to","bond_type"))%>%
    as.data.frame()%>%
    dplyr::mutate(
      bond_id = 1:n(),
      x = atom.data$x[from],
      xend = atom.data$x[to],
      y = atom.data$y[from],
      yend = atom.data$y[to]
    )%>%
    dplyr::mutate(xl = (xend-x),
                  yl = (yend - y),
                  x = x + bond.length.short*xl,
                  xend = xend - bond.length.short*xl,
                  y = y+bond.length.short*yl,
                  yend = yend - bond.length.short*yl)
  for (i in 1:nrow(bond.data)) {
    bond.data <- dplyr_copy_row(bond.data,
                                i,
                                bond.data$bond_type[i]-1)
  }

  lw <- 0.5*cex
  sw <- 0.5*cex
  col.bond <- "#666666"
  ggplot()+
    ### 3 bond
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw+2*sw+2*lw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = lw+2*sw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    ### 2 bond
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = 2*lw+sw)+
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = sw)+
    ### 1 bond
    geom_segment(data = filter(bond.data,bond_type == 1),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    geom_text(aes(x = median(range(atom.data$x)),
                  y = max(atom.data$y)+diff(range(atom.data$y))*0.5,
                  label = paste(sdf.formula,"\n",sdf.mz)),
              size = 2)+
    ylim(c(min(atom.data$y),max(atom.data$y)+diff(range(atom.data$y))*0.8))+
    xlim(expand_range(range(atom.data$x),multi = 0.2))+
    theme_void()->p


  if (show_ele) {
    p <- p+geom_text(data = atom.data,
                aes(x = x, y = y ,label = element),
                size = 2 *cex)
  }else{
    p <- p+geom_point(data = atom.data,
                     aes(x = x, y = y ),
                     size = 0.5*cex)
  }
  p
  return(p)
}

check_sdf <- function(sdf){

  atom.map.matrix <- atomcountMA(sdf)
  atom.map.matrix <- atom.map.matrix[,setdiff(colnames(atom.map.matrix),"0"),drop =F]
  id.atom <- apply(atom.map.matrix,1,sum)>1

}

check_smile <- function(smile){

  smile.sdf <- get_smile_sdf(smile)
  check_sdf(smile.sdf)

}

#' @export
get_sdf_formula <- function(sdfs){

  if (class(sdfs)=="SDF"    ) {
    sdfs <- ChemmineR::SDFset(list(sdfs))
  }
  sdfs.checked <- check_sdf(sdfs)
  sdfs.formula <- character()
  sdfs.formula[sdfs.checked] <- MF2(sdfs[sdfs.checked],addH=T)
  sdfs.formula <- MSCC::chemform_formate(sdfs.formula)
  return(sdfs.formula)
}

#' get_smiles_sdf
#'
#' @title Convert SMILES strings to SDF format
#' @description Converts one or more SMILES strings to SDF (Structure Data Format)
#'   objects using the ChemmineR package. Optionally canonicalizes the structures.
#'   Uses a precomputed mapping table to replace known SMILES with stored SDFs.
#' @param smiles Character vector of SMILES strings to convert.
#' @param smiles.id Optional character vector of IDs to assign to the resulting SDF objects.
#'   If NULL and `smiles` has names, those names are used; otherwise IDs are generated as "CMP001", etc.
#' @param canonicalize Logical indicating whether to canonicalize the SDF structures (default: TRUE).
#' @return An SDFset object (list of SDF objects) containing the molecular structures.
#' @export
get_smiles_sdf <- function(smiles,
                           smiles.id = names(smiles),
                           canonicalize = T){

  if (is.null(names(smiles))) {
    if (is.null(smiles.id)) {
      names(smiles) <- paste0("CMP",num2str(seq_along(smiles)))
    }else{
      names(smiles) <- smiles.id
    }

  }
  smiles.sdf <- suppressWarnings(
  ChemmineR::smiles2sdf(smiles)
  )

 #sdfs <- plyr::llply(smiles,.progress = "text",.fun = function(x){
 #  if (is.na(x)) return(NA)
 #  #print(x)
 #  ChemmineR::smiles2sdf(x)[[1]]
 #})
 #smiles.sdf <- ChemmineR::SDFset(sdfs)

  data(smiles_map)
  for (id in ChemmineR::cid(smiles_map)) {
    which(smiles==id)
    suppressWarnings(smiles.sdf[smiles==id] <- smiles_map[[id]])

  }


  if (canonicalize) {

    #sdfs <- plyr::llply(a,.progress = "text",.fun = function(x){
    #  if (is.na(x)) return(NA)
    #  #print(x)
    #  ChemmineR::canonicalize(x)
    #})
    #smiles.sdf <- ChemmineR::SDFset(sdfs)
    smiles.sdf <- ChemmineR::canonicalize(smiles.sdf)
  }
  #cid(smiles.sdf) <- smiles.id
  return(smiles.sdf)

}

#' @export
get_sdf_smiles <- function(sdf){
  if (inherits(sdf,"SDF")) {
    sdf <- ChemmineR::SDFset(list(sdf))
  }
  ChemmineR::sdf2smiles(sdf)

}

get_smile_formula <- function(smile){

  smile.sdf <- get_smiles_sdf(smile)
  smile.formula <- get_sdf_formula(smile.sdf)
  smile.formula <- case_when(smile=="O"~"H2O1",
                 smile=="[HH]"~"H2",
                 T~smile.formula)

  return(smile.formula)
}


vis_sdf_igraph_old <- function(sdf.igraph ,
                               show_id = F,
                               prob.border = NULL,
                               prob.fill = NULL,
                               highlight = NULL){

  sdf.igraph.temp <-sdf.igraph

  ### map prob to color and hight
  {

    ele <-  get_sdf_igraph_atom(sdf.igraph)
    #prob.fill <- prob.border <- runif(10,0,1)%>%`names<-`(sample(get_sdf_igraph_atom(sdf.igraph),10))
    if (is.numeric(highlight)|is.logical(highlight)) highlight <- ele[highlight]
    prob.border[highlight] <- 1
    col.border <- .get_vis_col(sdf.igraph,prob.border,
                               colramp(breaks = c(0,Inf,1),
                                       colors = c("#aaaaaa","#97C2FC","#2B7CE9")))
    col.fill <- .get_vis_col(sdf.igraph,prob.fill,
                             na.col = "#DDDDDD",
                             colramp(breaks = c(0,Inf,1),
                                     colors = c("#FFFFFF","#F7844F","#B20C26")))
    ele <- get_sdf_igraph_atom(sdf.igraph)

  }

  vda <- vdata(sdf.igraph.temp)<- vdata(sdf.igraph.temp)%>%
    dplyr::mutate(label = case_when(show_id~id,
                                    T~paste0(" ",atom," ")),
                  label = str_format_len(label),
                  font.size = case_when(show_id~20,T~40),
                  # font.multi= T,
                  # font.bold = T,
                  # font.bold.mod = "bold",
                  # font.bold.size = 500,
                  font.vadjust = 5,
                  # font.strokeWidth = 2,
                  #  font.strokeColor = "black",
                  font.align = "left",
                  borderWidth = 3,
                  color.background = col.fill[name],
                  color.border = col.border[name],
                  shape = "circle"
    )

  sdf.igraph.temp%>%
    visIgraph(idToLabel = F,
              type = "square")%>%
    visEdges(arrows = list(to = F),
             length = 2)

}




#' @export
get_sdf_igraph_atom <- function(ig,ele = "all"){

  vdf <- vdata(ig)
  if (ele== "all") {
    return(vdf$name)
  }else{
    vdf <- vdf %>%
      dplyr::filter(atom %in% ele)
    return(vdf$name)
  }

}



.get_highlight <- function(sdf.igraph,highlight){

}


vis_sdf <- function(sdf,show_id = F,...){

  sdf.igraph <- get_sdf_igraph(sdf)
  vis_sdf_igraph(sdf.igraph,show_id = show_id,...)


}




### to be removed
get_atom_id_from_parent <- function(parent.sdf.graph,
                                    product.sdf.graph){

  #parent.sdf.graph <- fragment.igraph[[1]]
  #product.sdf.graph <- fragment.igraph[[3]]

  old.root.id <- V(product.sdf.graph)$root_atom_id
  V(product.sdf.graph)$root_atom_id <- "unknown"
  ig <- intersection(parent.sdf.graph,
                     product.sdf.graph,
                     byname = F, keep.all.vertices = F)

  ig <- ig - V(ig)[atom_1!=atom_2]

  ig.vd <- vdata(ig)
  rownames(ig.vd) <- (ig.vd$name_2)

  #vis_sdf_igraph(parent.sdf.graph,show.label = F)
  #vis_sdf_igraph(product.sdf.graph,show.label = F)
  new.id <- ig.vd[V(product.sdf.graph)$name,]$root_atom_id_1
  if (!is.null(old.root.id)) {
    new.id <- case_when(is.na(new.id)~old.root.id,
              new.id == "unknown"~old.root.id,
              T~new.id)
  }

  V(product.sdf.graph)$root_atom_id <- new.id

  return(product.sdf.graph)

}


#' @title Map atoms between parent and product molecules
#' @description Performs atom mapping between two molecular structures using maximum common substructure (MCS)
#'   analysis. Computes probability mappings of atoms from parent to product, accounting for ring differences
#'   and bond type similarities.
#' @param sdf.parent An SDF object representing the parent molecule.
#' @param sdf.product An SDF object representing the product molecule.
#' @param ig.parent An igraph object representing the parent molecule's molecular graph. If NULL, it is computed from `sdf.parent`.
#' @param ig.product An igraph object representing the product molecule's molecular graph. If NULL, it is computed from `sdf.product`.
#' @param iso_ele Character string specifying the isotope element to consider for mapping (default: `"[13]C"`).
#' @param return.type Character string indicating the type of mapping to return. Either "most_prob" (default) returns a vector of most likely atom mappings, or "prob_matrix" returns a probability matrix.
#' @return For `return.type = "most_prob"`: a vector where names are atoms of the product and values are probabilities of mapping to atoms of the parent.
#'   For `return.type = "prob_matrix"`: a matrix with rows as atoms of the parent and columns as atoms of the product, containing probabilities.
#'   Both return values have an attribute "bond.score" indicating the bond similarity score.
#' @export
#'
get_atom_map <- function(sdf.parent,
                         sdf.product,
                         ig.parent = get_sdf_igraph(sdf.parent),
                         ig.product = get_sdf_igraph(sdf.product),
                         iso_ele = "[13]C",
                         return.type = c("most_prob","prob_matrix")){
  return.type <- match.arg(return.type)
  mcs <- fmcsR::fmcs(sdf.parent,
                     sdf.product,bu = 10)
  mcs.map <- MSCC::get_mcs_atom_map(mcs)
  mcs.map <- mcs.map.filter.duplicate(mcs.map,target_ele = MSCC::get_ele_uniso(iso_ele))
  product_atoms <- MSCC::atom(sdf.product)
  atom.map.matrix <- matrix(nrow = length(mcs.map),
                        ncol = length(product_atoms),
                        dimnames = list(seq_along(mcs.map),
                                        product_atoms))
  ring.diff <- length(ChemmineR::rings(sdf.parent))- length(ChemmineR::rings(sdf.product))
  bond.score <- rep(0,length(mcs.map))
  for (j in seq_along(mcs.map)) {
    this.map <- mcs.map[[j]]
    this.mapv <-this.map$mc1.atom
    names(this.mapv) <- this.map$mc2.atom
    this.mapv <- this.mapv[product_atoms]
    names(this.mapv) <- product_atoms


    ### ring re-assign
    {
      ring.solved <- F
      if (ring.diff&sum(is.na(this.mapv))) {
        ring.atom <- unname(unlist(ChemmineR::rings(sdf.parent)))
        ring.atom.to.assign <- ring.atom[!ring.atom%in% this.mapv]
        ring.atom.to.assign <- unique(ring.atom.to.assign)
        adj <- sapply(ring.atom.to.assign,function(x){
          #x <- ring.atom.to.assign
          x.adj <- names(igraph::V(ig.parent))[igraph::distances(ig.parent,x)==1]
          x.adj <- x.adj[x.adj%in%this.mapv]
          y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
          y.candi <-apply(igraph::distances(ig.product,y.adj),1,function(z){
            zz <- names(z)[which(z==1)]
             zz[!zz%in% names(na.omit(this.mapv))&
                  str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
            })%>%unlist()%>%unique()
          #message("ring candi: ",length(y.candi))
          unname(y.candi[1])
          })

        adj <- na.omit(unlist(adj))
        this.mapv[adj] <- names(adj)
        ring.solved <- ifelse(length(adj),T,F)
      }

    }

    ### ring nearest
    {
      if (ring.diff&sum(is.na(this.mapv))&ring.solved){
#
        ring.nearest.to.assign <-apply(igraph::distances(ig.parent,ring.atom),
                        1,function(z){
          zz <- names(z)[which(z==1)]
          zz[!zz%in% (na.omit(this.mapv))]
        })%>%unlist()%>%unique()

        adj <- sapply(ring.nearest.to.assign,function(x){
          #x <- ring.nearest.to.assign
          x.adj <- names(igraph::V(ig.parent))[igraph::distances(ig.parent,x)==1]
          x.adj <- x.adj[x.adj%in%this.mapv & x.adj%in%ring.atom]
          y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
          y.candi <-apply(igraph::distances(ig.product,y.adj),1,function(z){
            zz <- names(z)[which(z==1)]
            zz[!zz%in% names(na.omit(this.mapv))&
                 str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
          })%>%unlist()%>%unique()
          #message("ring nearest candi: ",length(y.candi))
          unname(y.candi[1])
        })

        adj <- na.omit(unlist(adj))
        this.mapv[adj] <- names(adj)


      }

    }

    ### non-match nearest
    {

     if (sum(is.na(this.mapv))) {
       non.match.to.assign <- setdiff(MSCC::atom(sdf.parent),this.mapv)
       adj <- sapply(non.match.to.assign,function(x){
         #x <- non.match.to.assign
         x.adj <- names(igraph::V(ig.parent))[igraph::distances(ig.parent,x)==1]
         x.adj <- x.adj[x.adj%in%this.mapv]
         y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
         y.candi <-apply(igraph::distances(ig.product,y.adj),1,function(z){
           zz <- names(z)[which(z==1)]
           zz[!zz%in% names(na.omit(this.mapv))&
                str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
         })%>%unlist()%>%unique()
         #message("non match nearest candi: ",length(y.candi))
         unname(y.candi[1])
       })

       adj <- na.omit(unlist(adj))
       this.mapv[adj] <- names(adj)

     }

    }

    ### bond diff
    {
      temp.map <- na.omit(this.mapv)
      temp.map.t <- make_vector(names(temp.map),temp.map)
      ig.sub <- igraph_filter_distance(ig.parent,from = temp.map,dis = 1)

      m1 <- igraph::as_adjacency_matrix(ig.sub,
                   attr = "bond_type",sparse = F)
      m1 <- m1+t(m1)
      m2 <- igraph::as_adjacency_matrix(ig.product,
                   attr = "bond_type",sparse = F)[names(temp.map),names(temp.map)]
      m2 <- m2 + t(m2)
      m2 <- get_matrix_value_fill_with_NA(
        m2,temp.map.t[rownames(m1)],temp.map.t[colnames(m1)])
      m2[is.na(m2)] <- 0
      bond.score[j] <- 1- sum((m1-m2)!=0)/sum(m1!=0)


    }
    #p<-vis_sdf_igraph_compare(ig.parent,ig.product,temp.map,names(temp.map),show.label = T)

    atom.map.matrix[j,] <- this.mapv


  }


  ### select map
  {
    atom.ele <- vdata(ig.product)$atom
    iso.atom.map.matrix <- atom.map.matrix[,atom.ele==MSCC::get_ele_uniso(iso_ele ),drop = F]
    atom.count <- apply(iso.atom.map.matrix,1,function(x)sum(!is.na(x)))
    full.mapped <- apply(iso.atom.map.matrix,1,function(x)sum(is.na(x))==0)
    selected <- which((atom.count+bond.score)==max(atom.count+bond.score,na.rm = T))

  }
  if (return.type == "most_prob"){
    map <-  apply(atom.map.matrix[selected,,drop=F],2,function(x){
      x <- na.omit(x)
      xp <- table(x)/length(x)
      parent_atoms <- MSCC::atom(sdf.parent)
      xp <- xp[parent_atoms]
      names(xp) <- parent_atoms
      xp[is.na(xp)] <- 0
      return(xp)
    })

    attributes(map)$bond.score <- mean(bond.score[selected])
  }


  if (return.type== "prob_matrix"){
    #apply(atom.map.matrix[atom.count==max(atom.count),,drop=F],2,function(x){
    map <-  apply(atom.map.matrix[full.mapped,,drop=F],2,function(x){
      x <- na.omit(x)
      xp <- table(x)/length(x)
      parent_atoms <- MSCC::atom(sdf.parent)
      xp <- xp[parent_atoms]
      names(xp) <- parent_atoms
      xp[is.na(xp)] <- 0
      return(xp)
    })
    attributes(map)$bond.score <- mean(bond.score)

  }



  return(map)

}

mcs.map.filter.duplicate <- function(mcs.map,target_ele = "C"){


  mcs.list <- lapply(mcs.map,function(x){

    x <- x %>%
      dplyr::filter(
        grepl(target_ele,mc1.atom)
      )%>%
      dplyr::arrange(mc1.idx)

    make_vector(x$mc2.atom,name = x$mc1.atom)

  })

  return(mcs.map[!duplicated(mcs.list)])
}



#' @title Visualize a molecule from SMILES string
#' @description Creates an interactive visualization of a molecular structure from a SMILES string.
#'   Converts the SMILES to an SDF object, then to an igraph representation, and generates an HTML widget
#'   using visNetwork. Optionally displays the molecular formula.
#' @param smiles A single SMILES string representing the molecule to visualize.
#' @param show.formula Logical indicating whether to display the molecular formula in the plot (default: TRUE).
#' @param show_id Logical indicating whether to show atom IDs as labels (default: TRUE). If FALSE, atom symbols are shown.
#' @param highlight Optional character vector of atom IDs to highlight in the visualization.
#' @return An HTML widget object (visNetwork) that can be rendered in RStudio viewer or browser.
#' @export
vis_smiles <- function(smiles,
                       show.formula = T,
                       show_id =T,
                       highlight =NULL){

  smiles.sdf <- get_smiles_sdf(smiles)[[1]]
  smiles.igraph <- get_sdf_igraph(smiles.sdf)
  smiles.vis <- vis_sdf_igraph(smiles.igraph,
                               show_id = show_id,
                               highlight = highlight)

  if (show.formula) {
    smiles.vis$x$main<- list(text = unname(ChemmineR::MF(smiles.sdf,addH = T)),
                             style = "text-align:center")
  }
  smiles.vis

}



#' @export
get_isopattern_score <- function(formula,
                                 mzs,
                                 int_matrix,
                                 ppm = 10){


  if (!length(formula)) return(NULL)
  formula.f <- factor(formula)
  iso_patterns <- lapply(levels(formula.f),
                        MSCC::chemform_isotopes_pattern_enviPat )
  iso_pattern <- iso_patterns[[1]]
  ip.score  <- lapply(iso_patterns,
         function(iso_pattern){
           if (nrow(iso_pattern)<=1) return(NA)
           iso_patterng <-iso_pattern %>%
             dplyr::ungroup()%>%
             dplyr::mutate(groupMz(x =m.z, ppm=ppm,return.type = "d"))%>%
             dplyr::group_by(mz.center)%>%
             dplyr::mutate(abundance=sum(abundance))%>%
             dplyr::distinct(mz.center,abundance)%>%
             dplyr::ungroup()

           id <- match_mz(mz1 = iso_patterng$mz.center,
                             mz2 = mzs,
                             mz.ppm  = ppm)
           iso.valm <- int_matrix[id,,drop = F]
           iso.ratio <- t(t(iso.valm)/iso.valm[1,])*100
           apply(iso.ratio , 2, function(iso.ratio.x){
             x <- iso_patterng$abundance[-1]
             y <- iso.ratio.x[-1]
             if (!length(x)) return(NA)
             if (all(is.na(y)))  return(0)
             y[is.na(y)] <- 0
             id.na <- is.na(x)|is.na(y)
             x <- x[!id.na]
             y <- y[!id.na]
             sum(x*y)^2/(sum(x^2L) *
                           sum(y^2L))
             1/exp(weighted.mean((abs(x-y)/x),w = x))
           })


         })

  ip.score <- sapply(ip.score,mean)
  ip.score <- ip.score[as.numeric(formula.f)]
  return(ip.score)
}


make_element_table <- function(){

  data("elem_table",package = "lc8")

  element_table <- elem_table%>%
    dplyr::mutate(element = MSCC::get_ele_uniso(element) ,
                  isotope = gsub("([0-9]+)", "[\\1]", isotope))%>%
    dplyr::group_by(element)%>%
    dplyr::mutate(is.isotope = case_when(
      abundance==max(abundance) ~F,
      T~T
    ),
    symbol = case_when(
      is.isotope~isotope,
      T~element
    ))



}


#' @export
is.isotope <- function(atoms){

  data(element_table)
  m <- make_vector(element_table$is.isotope,element_table$symbol )[atoms]
  dim(m) <- dim(atoms)
  dimnames(m) <- dimnames(atoms)
  return(m)
}




#' @export
MF2 <- function (x, ...){

  if (class(x) == "SDF")
    x <- as(x, "SDFset")
  propma <- atomcountMA(x, ...)
  propma <- propma[c(1, seq(along = propma[, 1])), ,drop = F]
  hillorder <- colnames(propma)
  names(hillorder) <- hillorder
  hillorder <- na.omit(unique(hillorder[c("C", "H", sort(hillorder))]))
  propma <- propma[, hillorder,drop = F]
  propma[propma == 1] <- ""
  MF <- paste(colnames(propma), t(propma), sep = "")
  propma <- matrix(MF, nrow = length(propma[, 1]), ncol = length(propma[1,
  ]), dimnames = list(rownames(propma), colnames(propma)),
  byrow = TRUE)
  MF <- seq(along = propma[, 1])
  names(MF) <- rownames(propma)
  zeroma <- matrix(grepl("[\\*A-Za-z]0$", propma), nrow = length(propma[,
                                                                        1]), ncol = length(propma[1, ]), dimnames = list(rownames(propma),
                                                                                                                         colnames(propma)))
  propma[zeroma] <- ""
  for (i in seq(along = MF)) {
    MF[i] <- paste(propma[i, ], collapse = "")
  }
  return(MF[-1])
}



#' @export
format_isotopologue <- function(x,
                                format = c("number","M","M+","+")){
  format <- match.arg(format)
  if (format == "+") format <- "M+"
  format <- match.arg(format,c("number","M","M+"))
  x <- str_extract_num(x)
  if(!length(x)) return(NULL)
  switch(
    format,
    "number" = x,
    "M" = paste0("M",x),
    "M+" = paste0("M+",x)
  )



}


#' @export
make_isotopologues_col <- function(n=10){


  suppressWarnings(
    cols <- c(ggsci::pal_npg()(10),ggsci::pal_bmj()(10))%>%na.omit()
  )
  #scales::show_col(cols)
  setNames(cols[1:(n+1)],format_isotopologue(0:n,format = "+"))

}


get_formula_from_CN_mz <- function(
    mz = 810.133057,
    CN_formula = "C23N7",
    charge = 1,
    ppm = 10
    ){


  c.count <- MSCC::get_formula_ele_count(CN_formula,"C")
  n.count <- MSCC::get_formula_ele_count(CN_formula,"N")
  mzo <- mz - MSCC::chemform_mz(CN_formula)

  lc8::mz_formula(
    Accurate_mass = mz,
    charge = charge,
    ppm = ppm,
    C_range = c.count:c.count,
    N_range = n.count:n.count,
    H_range = 0:floor(min(mzo, c.count * 2 + 2 + n.count + 5)),
    O_range = 0:floor(mzo/16),
    Cl_range = 0:floor(min(mzo/35,3)),
    P_range = 0:floor(mzo/31),
    S_range = 0:floor(mzo/32),
    Na_range = 0:floor(min(mzo/23,1)),
    K_range = 0:floor(min(mzo/16,1)),
    F_range = 0:floor(min(mzo/19,5)),
    Br_range = 0:floor(min(mzo/79,2)),
    I_range = 0:0,
    Si_range = 0:0,
    B_range = 0:0,
    Ca_range = 0:0,
    Cu_range = 0:0,
    Ni_range = 0:0,
    N_rule = T,
    Elem_ratio_rule = F,
    db_min = 0,
    db_max = 99,
    metal_ion = 0:3
  )
 #list(
 #  "H" = 0:floor(min(mzo, c.count * 2 + 2 + n.count)),
 #  "O"= 0:floor(mzo/16),
 #  "P"= 0:floor(mzo/31),
 #  "S"= 0:floor(mzo/32),

 #  "K"= 0:floor(min(mzo/16,1)),
 #  "Na"= 0:floor(min(mzo/23,1)),


 #  "F"= 0:floor(min(mzo/19,5)),
 #  "Cl"= 0:floor(min(mzo/35,3)),
 #  "Br"= 0:floor(min(mzo/79,2))
 #)->ele.limits

#
 # ele.matrix <- expand.grid(ele.limits)
 # tmz <- MSCC:::chemform_matrix_mz(ele.matrix,charge = charge)
 # idm <- match_mz_foverlaps(mzo,tmz)
#
 # ele.matrix[idm$ion2,]%>%
 #   MSCC:::chemform_from_ele_matrix()
#
}




