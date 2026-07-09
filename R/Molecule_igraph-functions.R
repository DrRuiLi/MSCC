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




get_Molecule_igraph_from_MSIPAtomMap <- function(msipAtomMap, fragment.id = 1){

  get_Molecule_igraph_from_sdf(msipAtomMap@fragment_sdf[[fragment.id]])
}

get_Molecule_igraph_from_cfmd <- function(cfmd, fragment.id = 1){
  .Deprecated("get_Molecule_igraph_from_MSIPAtomMap")
  get_Molecule_igraph_from_MSIPAtomMap(msipAtomMap = cfmd, fragment.id = fragment.id)
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

get_Molecule_igraph_MS2 <- function(Molecule_igraph,cfmd){


  isotopomers <- Molecule_igraph@isotopomer
  FG.map <- cfmd@fragment_group_map

  #### prob
  {
    FG.map[FG.map<0.5] <- 0
    FG.map[FG.map>=0.5] <- 1
  }

  data(element_table)
  isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
  ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
  isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
  dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
  dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)

  isotopomers.mz.diff <- isotopomers.mz.diff[,colnames(FG.map)]
  fg.isotopomers.mz.diff <- isotopomers.mz.diff %*% t(FG.map)


  ms2.data <- cbind( cfmd@fragment_group,t(fg.isotopomers.mz.diff))%>%
    dplyr::filter(fragment_mz > 0)%>%
    tidyr::pivot_longer(rownames(isotopomers.eles),names_to = "isotopomer",values_to = "mzdiff")%>%
    dplyr::mutate( isotopomers[match(isotopomer,isotopomers$isotopomer),
                               c("isotopologue","label","abundance")],
                   mz = fragment_mz +mzdiff)%>%
    dplyr::group_by(mz)%>%
    dplyr::mutate(intensity = sum(abundance))



  return(ms2.data)


}

test_fun <- function(x){
  x
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



#' @description Molecule atom transfer by atom map.
#' @title Auto map atom structure
#' @describeIn Molecule_atom_transfer
#' return a matrix, column as atom of `to`, row as multiple map, value as atom of `from`
#'
#' @param mol.ig.from Molecule_igraph
#' @param mol.ig.to Molecule_igraph
#'
#' @returns `matirx`, column as atom of `to`, row as multiple map, value as atom of `from`
#' @export
#'
get_Molecule_atom_transfer_by_atom_map <- function(mol.ig.from,
                                              mol.ig.to,
                                              target_ele = "ANY"){


  if (target_ele=="ANY") {

    data(element_table)
    target_ele <- element_table$element
  }

  sdf.parent <- mol.ig.from@sdf
  sdf.product <- mol.ig.to@sdf
  ig.parent <- mol.ig.from@igraph
  ig.product <- mol.ig.to@igraph

  ### mcs map
  {
    mcs <- fmcsR::fmcs(sdf.parent,
                       sdf.product,bu = 10)
    mcs.map <- MSCC::get_mcs_atom_map(mcs)
    #mcs.map <- mcs.map.filter.duplicate(mcs.map,target_ele = get_ele_uniso(iso_ele))
  }


  ### complement mcs map
  {
    atom.map.matrix <- matrix(nrow = length(mcs.map),
                              ncol = length(atom(sdf.product,element = target_ele)),
                              dimnames = list(seq_along(mcs.map),
                                              atom(sdf.product,element = target_ele)))
    ring.diff <- length(rings(sdf.parent))- length(rings(sdf.product))
    bond.score <- rep(0,length(mcs.map))
    for (j in seq_along(mcs.map)) {
      this.map <- mcs.map[[j]]
      this.mapv <-this.map$mc1.atom
      names(this.mapv) <- this.map$mc2.atom
      this.mapv <- this.mapv[atom(sdf.product)]
      names(this.mapv) <- atom(sdf.product)


      ### ring re-assign
      {
        ring.solved <- F
        if (ring.diff&sum(is.na(this.mapv))) {
          ring.atom <- unname(unlist(rings(sdf.parent)))
          ring.atom.to.assign <- ring.atom[!ring.atom%in% this.mapv]
          ring.atom.to.assign <- unique(ring.atom.to.assign)
          adj <- sapply(ring.atom.to.assign,function(x){
            #x <- ring.atom.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
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
          ring.nearest.to.assign <-apply(distances(ig.parent,ring.atom),
                                         1,function(z){
                                           zz <- names(z)[which(z==1)]
                                           zz[!zz%in% (na.omit(this.mapv))]
                                         })%>%unlist()%>%unique()

          adj <- sapply(ring.nearest.to.assign,function(x){
            #x <- ring.nearest.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv & x.adj%in%ring.atom]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
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
          non.match.to.assign <- setdiff(atom(sdf.parent),this.mapv)
          adj <- sapply(non.match.to.assign,function(x){
            #x <- non.match.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
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

      atom.map.matrix[j,] <- this.mapv[atom(sdf.product,element = target_ele)]


    }
    atom.map.matrix <- unique(atom.map.matrix)
  }


  ### class
  {

    atom.map.sum <- data.frame(
      from.smiles = mol.ig.from@molecule_info$smiles,
      to.smiles = mol.ig.to@molecule_info$smiles,
      transfer = paste0("Atom_map",num2str(1:nrow(atom.map.matrix))),
      source = "Atom_map"
    )

    rownames(atom.map.matrix) <- atom.map.sum$transfer


    new("Molecule_atom_transfer",
        transfer_def = atom.map.sum,
        transfer_matrix = atom.map.matrix)



  }

}

get_Reaction_atom_transfer_by_atom_map <- function(mol.ig.from,
                                                   mol.ig.to,
                                                   target_ele = "ANY",
                                                   equation){


  if (target_ele=="ANY") {

    data(element_table)
    target_ele <- element_table$element
  }


  ### atom index
  {
    atom.from <- lapply(1:length(mol.ig.from),function(x){
      paste0(names(mol.ig.from)[x],"_",atom(mol.ig.from[[x]],element = target_ele))
    })%>%unlist()%>%unname

    atom.to<- lapply(1:length(mol.ig.to),function(x){
      paste0(names(mol.ig.to)[x],"_",atom(mol.ig.to[[x]],element = target_ele))
    })%>%unlist()%>%unname

  }


  ### atom map
  {
    mat.df <- expand.grid(from = names(mol.ig.from),
                          to = names(mol.ig.to))
    plyr::mlply(mat.df,function(from,to){
      get_Molecule_atom_transfer_by_atom_map(
        mol.ig.from = mol.ig.from[[from]],
        mol.ig.to = mol.ig.to[[to]],
        target_ele = target_ele
      )
    })->mats

    mats.map <- sapply(1:nrow(mat.df),function(i){
      x <- mats[[i]]@transfer_matrix
      colnames(x) <- paste0( mat.df$to[i],"_" , colnames(x))
      x <- matrix(ifelse(is.na(x), NA, paste0(mat.df$from[i],"_", x)),
                  nrow = nrow(x), dimnames = dimnames(x))
    })

    combs <- sapply(mats,function(x){1:length(x)})%>%do.call(expand.grid,.)
    lapply(1:nrow(combs),function(i){

      maps <- sapply(1:ncol(combs),function(j) na.omit(mats.map[[j]][combs[i,j],]))
      map <- do.call(bind_rows,maps)

      return(map)
    })->a
  }






}

#' @export
get_Reaction_atom_transfer_by_RXNmapper <- function(mol.ig.from,
                                                    mol.ig.to,
                                                    target_ele = "ANY",
                                                    equation){

  if (target_ele=="ANY") {
    data(element_table)
    target_ele <- element_table$element
  }


  ### RXNmapper
  {

    rxn.result <- MSCC::RXNMapper_map(from.smiles = sapply(mol.ig.from,function(x)x@molecule_info$smiles),
                  to.smiles = sapply(mol.ig.to,function(x)x@molecule_info$smiles))
    if (all(is.na((rxn.result)))) return(NA)
    rxn.atom.id <- MSCC::RXNMapper_mapped_rxn_parse(rxn.result)
  }


  ### make atom_transfer
  {
    from.df <- lapply( names(mol.ig.from), function(i){
      make_vector(x = rxn.atom.id$from[[i]],
                  name = paste0(i,"_",atom(mol.ig.from[[i]])))
      data.frame(
        from.compound.id = mol.ig.from[[i]]@molecule_info$id,
        from.molecule.id = i,
        from.atom.id = atom(mol.ig.from[[i]]),
        from.rxn.id = rxn.atom.id$from[[i]]
      )%>%
        dplyr::filter(element(mol.ig.from[[i]]) %in% target_ele)
    })%>%do.call(bind_rows,.)

    to.df <- lapply( names(mol.ig.to), function(i){

      data.frame(
        to.compound.id = mol.ig.to[[i]]@molecule_info$id,
        to.molecule.id = i,
        to.atom.id = atom(mol.ig.to[[i]]),
        to.rxn.id = rxn.atom.id$to[[i]]
      )%>%
        dplyr::filter(element(mol.ig.to[[i]]) %in% target_ele)
    })%>%do.call(bind_rows,.)

    idx <- intersect(from.df$from.rxn.id,to.df$to.rxn.id)%>%sort()
    atom_transfer <- cbind(from.df[match(idx,from.df$from.rxn.id),],
          to.df[match(idx,to.df$to.rxn.id),])
    rownames(atom_transfer) <- NULL
  }


  ### make Reaction_atom_transfer
  {

    rat <- Reaction_atom_transfer()
    rat@atom_transfer <- atom_transfer


  }

  return(rat)


}

canonicalize_smiles <- function(smiles){

  smiles.sdf <- MSCC::get_smiles_sdf(smiles,canonicalize = T)
  smiles.canonicalized <-  ChemmineR::sdf2smiles(smiles.sdf)%>%
    as.character()
  names(smiles.canonicalized) <- names(smiles)
  return(smiles.canonicalized)
}


vis_Molecule_atom_transfer <- function(mat,id = 1,show_id = F){

  ### mol.ig
  {
    smiles.from <- unique(mat@transfer_def$from.smiles)
    smiles.to<- unique(mat@transfer_def$to.smiles)

    mig.from <- get_Molecule_igraph_from_smiles(smiles.from)
    mig.to <- get_Molecule_igraph_from_smiles(smiles.to)

  }

  ### transfer
  {

    transfer <- mat@transfer_matrix[id,]
    transfer <- transfer[!is.na(transfer)]
  }

  ### highlight
  {
    mig.from <- mig.from%>%
      Molecule_igraph_vis_format()%>%
      sdf_igraph_add_background_color(value = make_vector(0.5,transfer))
    mig.to <- mig.to%>%
      Molecule_igraph_vis_format()%>%
      sdf_igraph_add_background_color(value = make_vector(0.5,names(transfer)))
  }


  ### merge graph
  {

    ig.merge <- sdf_igraph_merge(mig.from,mig.to)%>%
      sdf_igraph_show_id(show_id)
    edata(ig.merge)$connect_type <- "bond"

    for (atom in names(transfer)) {
      ig.merge <- igraph::add_edges(ig.merge,
                                    edges = c(paste0("F_",transfer[atom] ),
                                      paste0("T_",atom)) ,
                                    attr = list(id = paste0("F_",transfer[atom],"-","T_",atom ),
                                                smooth.enabled = T,
                                                smooth.type = "dynamic",
                                                connect_type = "atom_transfer",
                                                width = 8,
                                                #color = "rgba(88, 203, 202, 0.5)",
                                                color = "#54BFBF",
                                                arrows.to = T,
                                                smooth.roundness = 0.5,
                                                shadow = T,
                                                dashes = I(list(c(10,20))))
                                    )
    }

  }


  ### Vis


  vis_igraph(ig.merge)




}



Molecule_atom_transfer_remove_duplicate <- function(mat,
                                          target_ele = "ANY"){


  if (target_ele=="ANY") {
    data(element_table)
    target_ele <- element_table$element
  }



  ### filter target_ele
  {

    smiles.to<- unique(mat@transfer_def$to.smiles)
    mig.to <- get_Molecule_igraph_from_smiles(smiles.to)
    target_atom <- atom(mig.to,target_ele)
  }


  ### duplicate
  {
    target_atom <- intersect(target_atom,colnames(mat@transfer_matrix))
    is.diplicate <- duplicated(mat@transfer_matrix[,target_atom,drop = F])
    mat@transfer_def <- mat@transfer_def %>%
      dplyr::filter(!is.diplicate)
    mat@transfer_matrix <- mat@transfer_matrix [mat@transfer_def$transfer,,drop = F]
  }


  return(mat)

}

Molecule_atom_transfer_reverse <- function(mat){

  mat@transfer_def <- mat@transfer_def %>%
    dplyr::mutate(tmp = to.smiles,
                  to.smiles = from.smiles,
                  from.smiles = tmp)%>%
    dplyr::select(!tmp)


  m <- mat@transfer_matrix
  m.list <- list()
  for (i in 1:nrow(m)) {

    x <- m[i,]
    x <- na.omit(x)
    m.list[[i]] <- make_vector(x = names(x),name = x)

  }
  m <- do.call(bind_rows,m.list)%>%as.matrix()
  rownames(m) <-rownames( mat@transfer_matrix)
  mat@transfer_matrix <- m

  return(mat)
}


get_mat_from_visGetEdges <- function(mat,id,visGetEdges){

  edge.data <- lapply(visGetEdges,`[`,c("from","to","id","connect_type"))%>%
    rbindlist(fill = T)
  edge.data <- edge.data%>%
    dplyr::filter(connect_type%in% c("atom_transfer","custom_atom_transfer"),
                  str_extract(from,"[FT]") =="F",
                  str_extract(to,"[FT]") =="T")

  atom.trans <- make_vector(x = sub("^F_", "", edge.data$from),
                            name =  sub("^T_", "", edge.data$to))
  atoms <- union(names(atom.trans),colnames(mat@transfer_matrix))
  x <- get_matrix_value_fill_with_NA(mat@transfer_matrix,colnames_vec = atoms)
  x[id,atoms] <- atom.trans[atoms]
  mat@transfer_matrix <- x
  mat <- Molecule_atom_transfer_remove_duplicate(mat)
  return(mat)
}


vis_Reaction_atom_transfer <- function(rat){


  ### mol.ig
  {
    from.migs <- MSCC::get_smiles_sdf(rat@reaction_info$from.smiles)%>%
      get_Molecule_igraph_from_sdf(id = names(rat@reaction_info$from.smiles))
    to.migs <- MSCC::get_smiles_sdf(rat@reaction_info$to.smiles)%>%
      get_Molecule_igraph_from_sdf(id = names(rat@reaction_info$to.smiles))

  }

  ### transfer
  {
    transfer <- rat@atom_transfer%>%
      dplyr::mutate(
        from.id = paste0(from.molecule.id,"_",from.atom.id),
        to.id = paste0(to.molecule.id,"_",to.atom.id)
      )
  }

  ### merge graph
  {
    from.mig <- sdf_igraph_merge_all(from.migs)
    to.mig <- sdf_igraph_merge_all(to.migs)

    ### merge from and to
    {
      xa <-vdata(from.mig)$x
      xb <- vdata(to.mig)$x

      vdata(from.mig)$x <- xa-max(xa)-diff(range(xa))/5
      vdata(to.mig)$x <- xb- min(xb)+diff(range(xb))/5

      nodes <- bind_rows(
        vdata(from.mig),
        vdata(to.mig)
      )
      #nodes$size <- 20
      edges <- rbind(
        edata(from.mig),
        edata(to.mig)
      )
      ig.merge <- igraph::graph_from_data_frame(edges,
                                                  vertices = nodes)%>%
        Molecule_igraph_vis_format%>%
        sdf_igraph_add_background_color(value = make_vector(0.5,
                                                            c(transfer$from.id,transfer$to.id)))
    }

    ### add transfer edge
    {
      for (i in 1:nrow(transfer)) {

        ig.merge <- igraph::add_edges(ig.merge,
                                      edges = c(transfer$from.id[i],
                                                transfer$to.id[i]) ,
                                      attr = list(id = paste0(
                                        transfer$from.id[i],"_",
                                        transfer$to.id[i]            ),
                                                  smooth.enabled = T,
                                                  smooth.type = "dynamic",
                                                  connect_type = "atom_transfer",
                                                  width = 8,
                                                  #color = "rgba(88, 203, 202, 0.5)",
                                                  color = "#54BFBF",
                                                  arrows.to = T,
                                                  smooth.roundness = 0.5,
                                                  shadow = T,
                                                  dashes = I(list(c(10,20))))
        )
      }

    }



  }




  ### Vis
  vis_igraph(ig.merge)



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

