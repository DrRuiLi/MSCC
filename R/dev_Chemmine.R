ggplot_sdf <- function(sdf,
                       cex = 1,
                       show_ele = F){
  stop("ggplot_sdf() no longer reads SDF objects. Use vis_smiles().", call. = FALSE)

  sdf.formula <- NULL
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
  stop("check_sdf() no longer reads SDF objects. Use check_smile().", call. = FALSE)
}

check_smile <- function(smile){

  smile.sdf <- get_smile_sdf(smile)
  check_sdf(smile.sdf)

}

#' @export
get_sdf_formula <- function(sdfs){
  if (!is.character(sdfs)) {
    stop("Pass SMILES to get_smile_formula().", call. = FALSE)
  }
  get_smile_formula(sdfs)
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
                           canonicalize = TRUE){
  stop(
    "get_smiles_sdf() is no longer an entry. Use rdkit_canonical_smiles() or get_Molecule_igraph_from_smiles().",
    call. = FALSE
  )
}

#' @export
get_sdf_smiles <- function(sdf){
  stop("get_sdf_smiles() is no longer an entry. Use rdkit_canonical_smiles().", call. = FALSE)
}

#' Molecular formula from SMILES
#'
#' @param smile Character SMILES string(s). Callers should unique the
#'   vector first when the same structure is repeated.
#' @return Character formula(s), formatted by [chemform_formate()].
#' @export
get_smile_formula <- function(smile){
  smile.formula <- vapply(as.character(smile), function(one) {
    if (is.na(one) || !nzchar(one)) {
      return(NA_character_)
    }
    MSCC::chemform_formate(rdkit_molecule_parts(one)$formula)
  }, character(1))
  smile.formula <- dplyr::case_when(
    smile == "O" ~ "H2O1",
    smile == "[HH]" ~ "H2",
    TRUE ~ unname(smile.formula)
  )
  smile.formula
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

  smiles.igraph <- get_sdf_igraph(smiles)
  smiles.vis <- vis_sdf_igraph(smiles.igraph,
                               show_id = show_id,
                               highlight = highlight)

  if (show.formula) {
    smiles.vis$x$main<- list(text = unname(get_smile_formula(smiles)),
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
  if (is.character(x)) {
    return(get_smile_formula(x))
  }
  stop("MF2 no longer reads SDF objects. Pass a SMILES string.", call. = FALSE)
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
