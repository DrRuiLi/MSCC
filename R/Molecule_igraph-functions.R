#' @export
get_Molecule_igraph_from_sdf <- function(sdf = sdfsample,
                                         id = paste0("CP",num2str(1:length(sdf)))) {


  if (class(sdf) == "SDF")
    Molecule_igraphs <- sdf_to_Molecule_igraph(sdf,id = id)


  if (class(sdf) == "SDFset") {
    Molecule_igraphs <- list()
    sdf.valid <- validSDF(sdf)
    for (i in 1:length(sdf)) {
      Molecule_igraphs[[i]] <- sdf_to_Molecule_igraph(sdf[[i]],id = id[i])
    }
    names(Molecule_igraphs) <- cid(sdf)
  }

  return(Molecule_igraphs)

}


sdf_to_Molecule_igraph <- function(sdf,
                                   id ) {


  bond.data <- bonds(sdf)
  bond.data <- bond.data[rownames(bond.data)!="0",,drop = F]
  atom.block <- atomblock(sdf) %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    dplyr::mutate(name = id,
                  bond.data,
                  element = atom,
                  x = C1,
                  y = C2,
                  .before = id)

  atom.ids <- make_vector(c(atom.block$id,"NA"),
                          c(as.character(seq_len(nrow(atom.block))),"0"
                            ))
  bond.block  <- bondblock(sdf) %>%
    as.data.frame() %>%
    dplyr::mutate(from = atom.ids[as.character(C1)],
                  to = atom.ids[as.character(C2)],
                  id = paste0(from,"_",to), ### identifier for visNet
                  bond_type = C3,
                  .before = C1)
  bond.block <- bond.block[!((bond.block$from=="NA")&
                               (bond.block$to=="NA")),]
  sdf.igraph <- igraph::graph_from_data_frame(bond.block, vertices = atom.block)

  isotopomer.df <-  make_vector(atom.block$element,atom.block$id)%>%
    as.list()%>%
    as.data.frame()%>%
    dplyr::mutate(isotopomer = paste0(rep(0,length(atom(sdf,element = "C"))),collapse = ""),
                  isotopologue = "M+0",
                  label = "base",
                  abundance = 1,
                  path = "",
                  .before = everything()
    )
  rownames(isotopomer.df) <- isotopomer.df$isotopomer
  Molecule_igraph <- new("Molecule_igraph",
                         sdf = sdf, igraph = sdf.igraph,isotopomer = isotopomer.df)
  Molecule_igraph@molecule_info$smiles <- unname(as.character(MSCC::get_sdf_smiles(sdf)))
  Molecule_igraph@molecule_info$id <- id
  return(Molecule_igraph)
}


#' @export
get_Molecule_igraph_from_smiles <- function(smiles = "NCC(O)=O",id ="A",canonicalize= T) {
  sdf <- MSCC::get_smiles_sdf(smiles,smiles.id = id,canonicalize = canonicalize)
  if (length(sdf)==1) sdf <- sdf[[1]]
  get_Molecule_igraph_from_sdf(sdf,id = id)
}




#' Add isotopomer to Molecule Igraph
#'
#' @description Adds an isotopomer to a Molecule_igraph object, updating its isotopomer data frame.
#'
#' @param Molecule_igraph A `Molecule_igraph` object.
#' @param isotopomer Isotopomer name. If NULL, a unique name is generated.
#' @param iso_vec Named vector mapping atom IDs to isotopic labels, e.g., c("C_1" = `"[13]C"`).
#' @param abundance Numeric abundance of the isotopomer.
#' @param path Path information (optional).
#' @param FSIS Fragment-specific isotopomer score (optional).
#'
#' @return The updated `Molecule_igraph` object.
#' @export
Molecule_igraph_add_isotopomer <- function(
    Molecule_igraph , isotopomer = NULL,iso_vec = NULL,abundance=NA,path = NA,FSIS = NA){

  if (identical(iso_vec,"all_C")) {
    iso_vec <- make_vector("[13]C",atom(glu.mi,"C"))
  }

  iso_label <- get_isotopomer_name(iso_vec)
  ele_vec <-make_vector(vdata(Molecule_igraph)$element,
                        atom(Molecule_igraph))
  ele_vec[names(iso_vec)] <- iso_vec
  isotopologue <- sum(MSCC::is.isotope(ele_vec),na.rm = T)
  isotopologue <- MSCC::format_isotopologue(isotopologue,"+")
  isotopomer.df <- Molecule_igraph@isotopomer
  if (is.null(isotopomer)){

    isotopomer <- setdiff(paste0(isotopologue,"_",seq_along(isotopomer.df$isotopomer)),
                         isotopomer.df$isotopomer)[1]
  }
  to.add <- data.frame(isotopomer = isotopomer,
                       isotopologue = isotopologue,
                       FSIS = FSIS,
                       label = iso_label ,
                       abundance = abundance,
                       path = path

                       )
  rownames(to.add) <- to.add$isotopomer
  isotopomer.df[isotopomer,names(to.add)] <- to.add
  isotopomer.df[isotopomer,atom(Molecule_igraph)] <- ele_vec

  isotopomer.df -> Molecule_igraph@isotopomer
  return(Molecule_igraph)


}


Molecule_igraph_add_isotopologue <- function(
    Molecule_igraph, isotopologue ,abundance = 0, target_ele = "[13]C" ,all_isotopomers = F){


  if (length(isotopologue)>1) {
    for (i in isotopologue) {
      Molecule_igraph <- Molecule_igraph_add_isotopologue(Molecule_igraph,
                                                          abundance = abundance,
                                                          isotopologue = i,
                                                          target_ele = target_ele)
    }
    return(Molecule_igraph)
  }


  mol.ig.atoms <- atom(Molecule_igraph,element = MSCC::get_ele_uniso(target_ele))
  atom.count <- length(mol.ig.atoms)
  iso.count <- MSCC::format_isotopologue(isotopologue,"n")
  if (all_isotopomers) {

    isom <- get_isotopomer_atom_matrix(atom.count,iso.count)

    for (i.isotopomer in 1:nrow(isom)) {
      Molecule_igraph <- Molecule_igraph_add_isotopomer(Molecule_igraph,
                                                        isotopomer = rownames(isom),
                                                        iso_vec = make_vector(target_ele,mol.ig.atoms[isom[i,]==1]),
                                                        abundance = abundance/nrow(isom),
                                                        FSIS = paste0(isotopologue,"_","FSIS"))
    }
  }else{

    Molecule_igraph <- Molecule_igraph_add_isotopomer(
      Molecule_igraph,
      isotopomer = paste0(isotopologue,"_","Fill"),
      iso_vec = make_vector(target_ele,sample( mol.ig.atoms, iso.count  )),
      abundance = abundance,
      FSIS = paste0(isotopologue,"_","FSIS"))

  }


  return(Molecule_igraph)

}


Molecule_igraph_fill_isotopologues <- function(Molecule_igraph,target_ele = "[13]C"){

  mol.ig.formula <- MSCC::get_sdf_formula(Molecule_igraph@sdf)
  istp.df <- Molecule_igraph@isotopomer

  atom.matrix <- as.matrix(istp.df[,atom(Molecule_igraph)])
  atom.count <- MSCC::get_formula_ele_count(mol.ig.formula,MSCC::get_ele_uniso(target_ele))

  iso.counts <- setdiff(0:atom.count,MSCC::format_isotopologue(istp.df$isotopologue))
  atoms <- atom(Molecule_igraph,MSCC::get_ele_uniso(target_ele))
  for (i.iso.count in iso.counts) {

    i.atom <- sample(atoms,i.iso.count)
    Molecule_igraph <- Molecule_igraph_add_isotopomer(
      Molecule_igraph,
      iso_vec = make_vector(target_ele,name = i.atom),
      abundance = 0,
      FSIS = MSCC::format_isotopologue(i.iso.count,"M")
      )

  }

  return(Molecule_igraph)
}


Molecule_igraph_filter_isotopomers <- function(Molecule_igraph,fraction_thresh = 0.001){

  istpm <- Molecule_igraph@isotopomer %>%
    dplyr::ungroup()%>%
    dplyr::mutate(abundance = abundance/sum(abundance,na.rm = T))%>%
    dplyr::filter(abundance >= fraction_thresh)

  Molecule_igraph@isotopomer <- Molecule_igraph@isotopomer%>%
    dplyr::filter(isotopomer %in% istpm$isotopomer)


  return(Molecule_igraph)
}


#' Label an isotopomer from its labeled atom ids.
#'
#' @param iso_vec Named character vector of isotope symbols.
#' @return A label string, or `"base"` when nothing is labeled.
#' @export
get_isotopomer_name <- function(iso_vec){

  iso_vec <- iso_vec[MSCC::is.isotope(iso_vec)]
  if (!length(iso_vec)) return("base")
  iso_label <- split(names(iso_vec),iso_vec)
  iso_label <- sapply(seq_along(iso_label),function(x){
    y <- iso_label[[x]]
    y <- sort(y)
    paste0("(",paste0(y,collapse = ","),")",names(iso_label)[x])
  })%>%paste0(collapse = ";")

}

Molecule_igraph_remove_isotopomer <-function(
  Molecule_igraph , isotopomer = NULL){

  x <- isotopomer
  Molecule_igraph@isotopomer <- Molecule_igraph@isotopomer%>%
    dplyr::filter(!isotopomer %in% x)

  return(Molecule_igraph)
}


get_Molecule_igraph_MS1 <- function(Molecule_igraph,polarity=1,adduct = NULL){

  if (is.null(adduct)  ) {
    adduct <- ifelse(polarity==0,"[M-H]-","[M+H]+")
  }


  ### ms spectra
  {

    data(element_table)
    mz.m0 <- MSCC::chemform_adduct(formula(Molecule_igraph),adduct = adduct)
    isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
    ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
    isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
    dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
    dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)
    mz.diff <- apply(isotopomers.mz.diff,1,sum)

    ms1.data <- Molecule_igraph@isotopomer%>%
      dplyr::mutate(mz = mz.m0+mz.diff)%>%
      dplyr::group_by(isotopologue)%>%
      dplyr::mutate(intensity = sum(abundance) )%>%
      dplyr::ungroup()%>%
      dplyr::select(any_of(c("isotopomer","isotopologue","label","abundance","mz","intensity")))

  nor.to <- ms1.data$intensity[which.min(ms1.data$mz)]
  ms1.data$intensity <-  ms1.data$intensity/nor.to*100
  }

  return(ms1.data)



}

Molecule_igraph_get_C_order <- function(Molecule_igraph){

  message("this function uncompleted")

  dism <- igraph::distances(Molecule_igraph@igraph)
  ele <- element(Molecule_igraph)
  ele.dis.R <- apply(dism,1,function(x){
    sum(x[ele!="C"])
  })
  ele.dis.C <- apply(dism,1,function(x){
    sum(x[ele=="C"])
  })
  ele.dis.R[ele=="C"]
  ele.dis.C[ele=="C"]

}


vis_Molecule_igraph <- function(Molecule_igraph,show_id = F){

  Molecule_igraph.formated <- Molecule_igraph%>%
    Molecule_igraph_vis_format()%>%
    sdf_igraph_show_id(show_id)



  vis_igraph(Molecule_igraph.formated) %>%
    visNetwork::visPhysics(enabled = F) %>%
    visNetwork::visOptions(width = "100%", height = "100%")

}


vis_Molecule_igraph_smiles <- function(smiles){
  mig <- smiles%>%
    get_Molecule_igraph_from_smiles()
    vis_Molecule_igraph(mig)
}


vis_Molecule_igraph_isotopomer <- function(Molecule_igraph,show_id = F,isotopomer = 1){

  Molecule_igraph.formated <- Molecule_igraph%>%
    Molecule_igraph_vis_format()%>%
    sdf_igraph_show_id(show_id)


  isotopomer.atoms <- Molecule_igraph@isotopomer[isotopomer,atom(Molecule_igraph)]%>%
    unlist()
  isotopomer.atoms <- isotopomer.atoms[MSCC::is.isotope(isotopomer.atoms)]
  labeled.atoms <-  make_vector(x = rep(1,length(MSCC::is.isotope(isotopomer.atoms))),
                                names(isotopomer.atoms))


  Molecule_igraph.formated%>%
    sdf_igraph_add_background_color(value = labeled.atoms)%>%
    vis_igraph() %>%
    visNetwork::visPhysics(enabled = F) %>%
    visNetwork::visOptions(width = "100%", height = "100%")

}

#' Format a molecule graph for visNetwork.
#'
#' @param Molecule_igraph A `Molecule_igraph`.
#' @return The same object with visNetwork vertex and edge fields.
#' @export
Molecule_igraph_vis_format <- function(Molecule_igraph){

  vdata(Molecule_igraph) <-  vdata(Molecule_igraph)%>%
   # dplyr::group_by(atom)%>%
    dplyr::mutate(
      label = paste0(" ",atom," "),
      x =  ( x-mean(x))*100,
      y = ( y-mean(y)) *100,
      font.size = 30,
      font.vadjust = 5,
      font.align = "center",
      color.border = "#AAAAAA",
      color.background = "#FFFFFF",
      borderWidth = 2,
      shape = "circle",
      physics = F)%>%
    dplyr::ungroup()


  edata(Molecule_igraph)  <- edata(Molecule_igraph)%>%
    dplyr::mutate(smooth = FALSE,
                  #arrows = FALSE,
                  width = 10*bond_type)

  return(Molecule_igraph)

}


canonicalize_smiles <- function(smiles){

  smiles.sdf <- MSCC::get_smiles_sdf(smiles,canonicalize = T)
  smiles.canonicalized <-  ChemmineR::sdf2smiles(smiles.sdf)%>%
    as.character()
  names(smiles.canonicalized) <- names(smiles)
  return(smiles.canonicalized)
}



#' @export
is_labeled <- function(mol.ig,target_ele = "C"){

  isotopomers.data <- mol.ig@isotopomer
  isotopomers.matrix <- isotopomers.data[,atom(mol.ig,target_ele)]%>%
    as.matrix()


  labeled.matrix <- MSCC::is.isotope(isotopomers.matrix)
  any(apply(labeled.matrix,1,sum,na.rm = T)>0)


}


plot_Molecule_igraph <- function(mol.ig,show_id = F,size = 2){


  tpf <- tempfile(fileext = ".html")
  vis <- vis_Molecule_igraph(mol.ig,show_id = show_id)
  htmlwidgets::saveWidget(vis, tpf, selfcontained = TRUE)
  #open_file(tpf)

  tpfpng <- tempfile(fileext = ".png")
  webshot2::webshot(tpf, tpfpng,quiet = T)
  #open_file(tpfpng)

  ggplot_from_img(tpfpng,size = size)

}


plot_Molecular_igraphs_isotopomer_bar <- function(mol.igs){

  mol.igs.itpms.list <- lapply(mol.igs,function(x){
    x@isotopomer
  })
  mol.igs.itpms.df <- data.table::rbindlist(mol.igs.itpms.list,idcol = "sample",fill = T)%>%
    dplyr::group_by(sample)%>%
    dplyr::mutate(abundance = abundance/sum(abundance))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(isotopologue , isotopomer)%>%
    dplyr::mutate(isotopomer = factor(isotopomer,levels = unique(isotopomer))  )%>%
    tidyr::complete(sample,isotopomer)
  mol.igs.itpms.df$abundance[is.na(mol.igs.itpms.df$abundance)] <- 0

  if (length(unique(mol.igs.itpms.df$isotopomer))<2) {
    return(NULL)
  }

  ggplot(mol.igs.itpms.df)+
    geom_bar(
      aes(x = isotopomer , y = abundance, fill = sample),
      stat = "identity",position = "dodge"
    )+
    ggsci::scale_fill_aaas()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -60,hjust = 0,vjust = 0))

}


plot_Molecular_igraph_isotopomer_circle <- function(mol.ig){




  isotopomer.data <- mol.ig@isotopomer
  if (nrow(isotopomer.data)==1) {
    isotopomer.data$FSIS <- "M0"
  }
  isotopomer.data <-isotopomer.data%>%
    dplyr::select(isotopomer,isotopologue,abundance,FSIS)%>%
    dplyr::mutate(abundance = abundance/sum(abundance),
                  r = rank(-abundance, ties.method = "random"),
                  others = F
                  #others = abundance < thresh

                  )%>%
    dplyr::group_by(others)%>%
    dplyr::mutate(
      abundance = case_when(
        others~sum(abundance),
        T~abundance
      ),
      isotopomer = case_when(
        others~"Others",
        #grepl("Fill",isotopomer)~"Unknown",
        T~ isotopomer
      ),
      isotopologue = case_when(
        others~"Others",
        T~ isotopologue
      ),
    )%>%
    dplyr::ungroup()%>%
    #dplyr::filter(r <= topN+1)%>%
    dplyr::filter(!others)%>%
    dplyr::arrange(isotopologue,abundance)%>%
    dplyr::mutate(isotopomer = factor(isotopomer,levels = (unique(isotopomer))))




  split.col <- function(col,idx){

    if(length(idx)==1) return(col)
    col <- col[1]
    col_fun <- circlize::colorRamp2(breaks = c(-1,max(idx)),
                                    colors = c("white",col))
    col_fun(idx)
  }

  all.istl <- unique(na.omit(isotopomer.data$isotopologue))
  col.isotopologues <- MSCC::make_isotopologues_col(max(MSCC::format_isotopologue(isotopomer.data$isotopologue,"n"),na.rm = T))
  #col.isotopologues["Others"] <- "#999999"

  isotopomer.data <- isotopomer.data%>%
    dplyr::group_by(isotopologue)%>%
    dplyr::mutate(col = col.isotopologues[isotopologue],
                  col = split.col(col,1:n()))%>%
    dplyr::ungroup()

  fsis.data <- isotopomer.data%>%
    dplyr::arrange(rev(isotopomer) )%>%
    dplyr::group_by(FSIS)%>%
    dplyr::mutate(fsis.r = sum(abundance))%>%
    dplyr::ungroup()%>%
    dplyr::distinct(FSIS,fsis.r)%>%
    dplyr::mutate(x = 1,
                  y = cumsum(fsis.r),
                  ys = c(0,head(y,-1))
                  )

  isotpologues.data <- isotopomer.data%>%
    dplyr::group_by(isotopologue)%>%
    dplyr::summarise(abundance = sum(abundance))%>%
    dplyr::arrange(MSCC::format_isotopologue(isotopologue,"n"))%>%
    dplyr::mutate(isotopologue = factor(isotopologue,levels = (unique(isotopologue))))

  p <- ggplot()+
    geom_bar(aes( x= 1.7, y = abundance,
                  fill = isotopologue,group = isotopologue),
             col = "black",
             data = isotpologues.data,
             linewidth = 1,width = 0.2,
             stat = "identity")+
    scale_fill_manual(
      name = "Isotopologues",
      values = col.isotopologues)+
    guides(fill = guide_legend(ncol = 1,order = 1)) +
    ggnewscale::new_scale_fill() +


    geom_bar(aes( x= 1, y = abundance,
                  #col = isotopologue,
                  fill = isotopomer),
             data = isotopomer.data,
             col = "white",
             linewidth = 0.2,width = 1,
             stat = "identity")+
    scale_fill_manual(
      name = "Isotopomers",
      values =setNames(isotopomer.data$col,isotopomer.data$isotopomer)
    )+
    geom_rect(data = fsis.data,
              aes(xmin =0.5,xmax = 1.5,
                  ymin = ys, ymax = y),
              fill = "transparent",col = "black",linewidth = 1)+

    #geom_vline(xintercept = c(0.525,1.475))+
    guides(fill = guide_legend(ncol = 2)) +
    xlim(c(0.3,2))+
    coord_polar(theta = "y",direction = -1)+
    theme_void(base_size  = 10)

  p





}

