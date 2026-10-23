#' @export
RXNMapper_map <- function(from.smiles ,
                           to.smiles,
                           RXNMapper = MSCC::get_RXNMapper()
                           ){

  #from.smiles.format <- sub(pattern = "\\.","-",from.smiles)
  #to.smiles.format  <- sub(pattern = "\\.","-",to.smiles)

  from.string <- paste0(from.smiles,collapse = ".")
  to.string <- paste0(to.smiles,collapse = ".")

  req <- paste0(from.string, ">>",to.string)
  rxns <- c(req,"C>>C")
  if (any(grepl("#",x = rxns))) {
    message(crayon::red("RXNMapper_map parsing with ###, remove it"))
    rxns <- gsub(pattern = "#","",rxns)
    from.smiles <- gsub(pattern = "#","",from.smiles)
    to.smiles <- gsub(pattern = "#","",to.smiles)

  }
  if (nchar(req)>511) {
    message(crayon::red("Reaction SMILES has too many tokens, should be at most 512."))
    return(NA)
  }

  rxn.result <- RXNMapper(rxns,detailed_output = T)[[1]]
  rxn.result$from <- from.smiles
  rxn.result$to <- to.smiles

  return(rxn.result)

}


#' @export
RXNMapper_mapped_rxn_parse <- function(rxn.result){


  rxn.result.split <- str_split(rxn.result$mapped_rxn,">>")[[1]]%>%
    sapply(function(x){
      str_split(x,"\\.")
    })
  names(rxn.result.split ) <- c("from","to")

  ### smiles combine, some of smiles contain . will be regarded as multiple compound
  {
    if (length(rxn.result.split$from) != length(rxn.result$from)) {
      cp.count <- str_count(rxn.result$from,"\\.")+1
      in.formula <- get_smile_formula(rxn.result$from)
      out.formula <- get_smile_formula(rxn.result.split$from)

      a <- which(out.formula %in%  in.formula[which(cp.count==1)])
      for (i in which(cp.count>1)) {

        b <- in.formula[i]
        comb <- combn(setdiff(seq_along(out.formula),a),cp.count[i])
        for (j in 1:ncol(comb)) {

          c <- do.call(MSCC::chemform_calc,
                  c(as.list(out.formula[comb[,j]]),return= "chemform"))
          if (sum(MSCC::chemform_calc(in.formula[i],c,"-"))==0) {
            d <- rxn.result.split$from[comb[,j]]
            rxn.result.split$from[comb[,j]] <- NA
            rxn.result.split$from <- c(rxn.result.split$from,
                                       paste0(d,collapse =  "."))
            rxn.result.split$from  <- na.omit(rxn.result.split$from )
          }
        }

      }

    }

    if (length(rxn.result.split$to) != length(rxn.result$to)) {
      cp.count <- str_count(rxn.result$to,"\\.")+1
      in.formula <- get_smile_formula(rxn.result$to)
      out.formula <- get_smile_formula(rxn.result.split$to)

      a <- which(out.formula %in%  in.formula[which(cp.count==1)])
      for (i in which(cp.count>1)) {

        b <- in.formula[i]
        comb <- combn(setdiff(seq_along(out.formula),a),cp.count[i])
        for (j in 1:ncol(comb)) {

          c <- do.call(MSCC::chemform_calc,
                       c(as.list(out.formula[comb[,j]]),return= "chemform"))
          if (sum(MSCC::chemform_calc(in.formula[i],c,"-"))==0) {
            d <- rxn.result.split$to[comb[,j]]
            rxn.result.split$to[comb[,j]] <- NA
            rxn.result.split$to <- c(rxn.result.split$to,
                                       paste0(d,collapse =  "."))
            rxn.result.split$to  <- na.omit(rxn.result.split$to )
          }
        }

      }

    }
  }


  rxn.canonical.id <- lapply(rxn.result.split,function(x){
    lapply(x, function(y){
      parsed <- rdkit_rxn_ids(y)
      ids <- parsed$rxn_id
      stats::setNames(
        vapply(ids, function(one) {
          if (is.null(one) || !length(one)) NA_integer_ else as.integer(one[[1]])
        }, integer(1)),
        names(ids)
      )
    })
  })

  ### asssign id
  {
    for(i in c("from","to")){

      x <- canonicalize_smiles(rxn.result.split[[i]])
      y <- canonicalize_smiles(rxn.result[[i]])
      idx <- match(x,y)
      while(anyDuplicated(idx)){
        id <- which(duplicated(idx))[1]
        y[idx[id]] <- NA
        idx[id] <-match(x[id],y)
      }
      names(rxn.canonical.id[[i]]) <- names(rxn.result[[i]])[idx]
      rxn.canonical.id[[i]] <- rxn.canonical.id[[i]][names(rxn.result[[i]])]

    }


  }

  return(rxn.canonical.id)


}



