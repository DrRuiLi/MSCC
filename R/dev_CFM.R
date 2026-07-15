CFM <- function(docker_run_param = "--rm",
                  result_dir = tempdir(),
                  cfmid_cmd = "cfm-id"
){
  paste0("docker run ",docker_run_param," ",
         " -v ",result_dir,":/cfmid",
         "  wishartlab/cfmid ",
         cfmid_cmd)->cmd
  return(cmd)

}


#' Normalize CFM adduct input to \code{"[M+H]+"} / \code{"[M-H]-"}
#'
#' Accepts polarity shorthand: \code{0}/\code{"0"}/\code{"-"} -> \code{"[M-H]-"};
#' \code{1}/\code{"1"}/\code{"+"} -> \code{"[M+H]+"}.
#' @keywords internal
.normalize_cfm_param_adduct <- function(param_adduct) {
  if (length(param_adduct) != 1L) {
    stop("param_adduct must be length 1", call. = FALSE)
  }
  if (is.numeric(param_adduct) || is.integer(param_adduct)) {
    if (param_adduct == 0) return("[M-H]-")
    if (param_adduct == 1) return("[M+H]+")
    stop("param_adduct numeric must be 0 (negative) or 1 (positive)", call. = FALSE)
  }
  ad <- as.character(param_adduct)[[1]]
  if (identical(ad, "-") || identical(ad, "0")) return("[M-H]-")
  if (identical(ad, "+") || identical(ad, "1")) return("[M+H]+")
  ad
}


CFM_get_param_config <- function(adduct = c("[M+H]+","[M-H]-"),
                                   param = T,
                                   config = T){
  adduct <- .normalize_cfm_param_adduct(adduct)
  adduct = match.arg(adduct, choices = c("[M+H]+", "[M-H]-"))

  param_file <- switch (adduct,
                        "[M+H]+" = " /trained_models_cfmid4.0/[M+H]+/param_output.log ",
                        "[M-H]-" =  " /trained_models_cfmid4.0/[M-H]-/param_output.log ",
                        " none "
  )
  config_file <- switch (adduct,
                        "[M+H]+" = " /trained_models_cfmid4.0/[M+H]+/param_config.txt ",
                        "[M-H]-" = " /trained_models_cfmid4.0/[M-H]-/param_config.txt ",
                        " none "
  )
  param_file <- ifelse(param,param_file," none ")
  config_file <- ifelse(config,config_file," none ")
  paste0(param_file,config_file)
}







setClass("CFM_data",
         slots = list(
           peak_assignment = "data.frame",
           fragment_define = "data.frame",
           fragment_transition = "data.frame"
         ),
         prototype = list(
           peak_assignment = data.frame(),
           fragment_define = data.frame(),
           fragment_transition = data.frame()
         ))


setMethod("show",signature ="CFM_data",definition =
            function(object){
              polarity_str <- if(nrow(object@fragment_define) > 0 && "polarity" %in% colnames(object@fragment_define)) {
                pol <- object@fragment_define$polarity[1]
                if(!is.na(pol)) {
                  if(pol == 0) "negative" else "positive"
                } else {
                  "unknown"
                }
              } else {
                "unknown"
              }
              message("CFM_data (", polarity_str, ") with ",nrow(object@fragment_define)," fragment")
            } )


#' Get Polarity Suffix
#' @title Get Polarity Suffix
#' @description Returns the suffix to append to fragment IDs based on polarity.
#' @param polarity Numeric polarity value (0 for negative, 1 for positive)
#' @return Character suffix "_0" for negative, "_1" for positive
#' @export
get_polarity_suffix <- function(polarity) {
  if (is.na(polarity)) return("")
  return(ifelse(polarity == 0, "_0", "_1"))
}


#' Predict Mass Spectra using CFM-ID
#' @title Predict Mass Spectra using CFM-ID
#' @description Predicts mass spectra fragments for a given molecule using the CFM-ID algorithm via Docker.
#' This function sends a SMILES or InChI string to the CFM-ID Docker container and returns predicted
#' fragment spectra at multiple collision energies.
#' @describeIn CFM predict
#' @param smiles_or_inchi_or_file SMILES string, InChI string, or path to a file containing molecular structure
#' @param prob_thresh Probability threshold for fragment prediction (default: 0.001)
#' @param param_adduct Adduct type for prediction, e.g., "\[M+H\]+" or "\[M-H\]-"
#'   (default: "\[M+H\]+"). Also accepts \code{0}/\code{1} (0 = "\[M-H\]-", 1 = "\[M+H\]+")
#'   and \code{"-"}/\code{"+"} as shorthand.
#' @param annotate_fragments Logical, whether to annotate fragments (default: 1)
#' @param output_file_or_dir Path to save results, or NULL to return results in memory (default: NULL)
#' @param apply_postproc Logical, whether to apply post-processing (default: 0)
#' @param suppress_exceptions Logical, whether to suppress exceptions (default: 1)
#'
#' @return A CFM_data object containing predicted spectra data, or path to output file if output_file_or_dir is specified
#' @export

CFM_predict <- function(smiles_or_inchi_or_file = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                          prob_thresh = 0.001,
                          param_adduct = "[M+H]+",
                          annotate_fragments = 1,
                          output_file_or_dir = NULL,
                          apply_postproc = 0,
                          suppress_exceptions = 1){

  param_adduct <- .normalize_cfm_param_adduct(param_adduct)

  out.file <- ifelse(is.null(output_file_or_dir),
                               tempfile(),
                               normalizePath(output_file_or_dir,mustWork = F))

  model_param = CFM_get_param_config(param_adduct)
  CFM(result_dir = dirname(out.file),
        cfmid_cmd = paste0(
          "cfm-predict ",
          smiles_or_inchi_or_file," ",
          prob_thresh," ",
          model_param,
          " 1 ",
          basename(out.file), " ",
          apply_postproc," ",
          suppress_exceptions

        ))->cmd

  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    polarity <- MSCC::get_polarity_from_adduct(param_adduct)
    cfm.result <- read_CFM_predict_result(out.file, polarity = polarity)
    return(invisible(cfm.result))
  }else{
    #message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }
}

#' Annotate Mass Spectra using CFM-ID
#' @title Annotate Mass Spectra using CFM-ID
#' @description Annotates experimental mass spectra with predicted fragments using the CFM-ID algorithm via Docker.
#' This function matches experimental spectra to predicted fragments based on mass tolerance.
#' @describeIn CFM annotate
#'
#' @param smiles_or_inchi SMILES string or InChI string of the molecule to annotate
#' @param spectrum_file Path to spectrum file (mgf, msp, etc.) or a Spectra object containing experimental spectra
#' @param id Identifier for the compound (default: "AN_ID")
#' @param ppm_mass_tol Mass tolerance in parts per million for matching (default: 5.0)
#' @param abs_mass_tol Absolute mass tolerance in m/z units (default: 0)
#' @param param_adduct Adduct type for annotation, e.g., "\[M+H\]+" or "\[M-H\]-"
#'   (default: "\[M+H\]+"). Also accepts \code{0}/\code{1} and \code{"-"}/\code{"+"}.
#' @param output_file Path to save annotation results, or NULL to return results in memory (default: NULL)
#' @param ... Additional arguments passed to underlying functions
#'
#' @return A CFM_data object containing annotated spectra data, or path to output file if output_file is specified
#'
#' @export
#'

CFM_annotate<- function(smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                          spectrum_file = NULL,
                          id = "AN_ID",
                          ppm_mass_tol = 5.0,
                          abs_mass_tol = 0,
                          param_adduct = "[M+H]+",
                          output_file = NULL,...){

  param_adduct <- .normalize_cfm_param_adduct(param_adduct)

  out.file <- ifelse(is.null(output_file),
                     tempfile(),
                     normalizePath(output_file,mustWork = F))
  win.dir<- dirname(out.file)
  if ("Spectra" %in% class(spectrum_file) ) {
    spectrum_file <- export_Spectra_peak_list_for_cfm(spectrum_file,tempfile())
  }
  file.copy(spectrum_file,
            paste0(win.dir,"/",basename(spectrum_file)))
  model_param = CFM_get_param_config(param_adduct,param = T,config = T)

  CFM(result_dir = win.dir,
        cfmid_cmd = paste0(
          "cfm-annotate ",
          smiles_or_inchi," ",
          basename(spectrum_file)," ",
          id," ",
          ppm_mass_tol, " ",
          abs_mass_tol," ",
          model_param,
          basename(out.file)
        ))->cmd
  cmd
  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file)) {
    polarity <- MSCC::get_polarity_from_adduct(param_adduct)
    return(invisible(read_CFM_annotate_result(out.file, polarity = polarity)))
  }else{
    message("Result export to: \n",crayon::red(output_file))
    return(invisible(out.file))
  }
}


#' Annotate Fragments using CFM Fraggen (Deprecated)
#' @title Annotate Fragments using CFM Fraggen
#' @describeIn CFM fraggen and annotate
#' @description
#' Deprecated. Use \code{\link{CFM_annotate_by_predict}} instead.
#' This function generates CFM_data by generating fragment graphs with CFM_fraggen,
#' not relying on Spectra input. It creates a synthetic peak assignment based on
#' the generated fragments.
#'
#' @note
#' Filters fragments with mz < 30 to avoid errors when converting SMILES to SDF.
#' This function is deprecated and will be removed in future versions.
#'
#' @param smiles_or_inchi SMILES string or InChI string of the molecule
#' @param spectrum_file Spectrum file path (currently unused, default: NULL)
#' @param max_depth Maximum depth for fragment generation (default: 1)
#' @param id Identifier for the compound (default: "AN_ID")
#' @param ppm_mass_tol Mass tolerance in ppm for matching (default: 5.0)
#' @param abs_mass_tol Absolute mass tolerance in m/z units (default: 0)
#' @param param_adduct Adduct type, e.g., "\[M+H\]+" or "\[M-H\]-" (default: "\[M+H\]+")
#' @param output_file Path to save results, or NULL to return results in memory (default: NULL)
#' @param ... Additional arguments
#'
#' @return A CFM_data object containing fragment definitions and peak assignments
#' @export
#'
CFM_annotate_by_fraggen <- function(
    smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
    spectrum_file = NULL,
    max_depth = 1,
    id = "AN_ID",
    ppm_mass_tol = 5.0,
    abs_mass_tol = 0,
    param_adduct = "[M+H]+",
    output_file = NULL,...){

  .Deprecated("CFM_annotate_by_predict")
  cfm.fragen <- CFM_fraggen(smiles_or_inchi,
                            max_depth = max_depth,
                            param_adduct = param_adduct)


  peak.data <- cfm.fragen@fragment_define%>%
    dplyr::filter(!intermediate)%>%
    dplyr::mutate(energy = "energy0",
                  intensity = NA,
                  fragment_score= NA)%>%
    dplyr::add_row(energy = c(rep("energy1",nrow(.)),
                              rep("energy2",nrow(.))),
                   fragment_mz=rep(.$fragment_mz,2),
                   fragment_id = rep(.$fragment_id,2))%>%
    dplyr::select(energy,
                  mz =fragment_mz,
                  intensity,
                  fragment_id,
                  fragment_score)%>%
    remove_rownames()

  cfm.fragen@peak_assignment <- peak.data
  return(cfm.fragen)


}


#' Convert CFM_data peak assignments to a Spectra object
#'
#' @title Convert CFM Data to Spectra Object
#' @description Converts CFM (Competitive Fragmentation Modeling) peak assignment
#'   data to a \code{Spectra} object with spectra at three collision energy levels
#'   (10, 20, 40 eV for energy0/1/2).
#' @param cfmd A \code{CFM_data} object containing \code{peak_assignment}.
#'   Typically from \code{\link{read_CFM_predict_result}} /
#'   \code{\link{CFM_annotate_by_predict}}.
#'
#' @return A \code{Spectra} object containing MS/MS spectra at three CE levels.
#' @export
get_CFM_data_Spectra <- function(cfmd) {
  if (!methods::is(cfmd, "CFM_data")) {
    stop("`cfmd` must be a CFM_data object.", call. = FALSE)
  }
  if (!requireNamespace("Spectra", quietly = TRUE)) {
    stop("Package Spectra is required to build CFM Spectra.", call. = FALSE)
  }
  cfm.df <- cfmd@peak_assignment
  if (is.null(cfm.df) || !nrow(cfm.df)) {
    return(Spectra::Spectra())
  }
  sp.list <- list()
  for (i in 0:2) {
    this.sp.data <- cfm.df %>%
      dplyr::filter(energy == paste0("energy", i)) %>%
      dplyr::distinct(mz, intensity)
    this.sp.df <- S4Vectors::DataFrame(
      collisionEnergy = switch(as.character(i),
                               "0" = 10,
                               "1" = 20,
                               "2" = 40)
    )
    this.sp.df$mz <- list(this.sp.data$mz)
    this.sp.df$intensity <- list(this.sp.data$intensity)
    sp.list[[paste0("energy", i)]] <- Spectra::Spectra(this.sp.df)
  }
  do.call(c, unname(sp.list))
}

#' Export a Spectra object to a CFM-ID peak-list text file
#'
#' Writes energy0/1/2 blocks (CE 10/20/40) for \code{\link{CFM_annotate}}.
#'
#' @param sp A \code{Spectra} object with \code{collisionEnergy}.
#' @param file Output file path.
#' @return Invisibly, \code{file}.
#' @export
export_Spectra_peak_list_for_cfm <- function(sp, file) {
  if (!requireNamespace("Spectra", quietly = TRUE)) {
    stop("Package Spectra is required.", call. = FALSE)
  }
  ce <- Spectra::collisionEnergy(sp)
  readr::write_lines(character(0), file = file)
  for (this_ce in unique(ce)) {
    sp.this <- sp[ce == this_ce]
    energy.type <- switch(
      as.character(this_ce),
      "10" = "energy0",
      "20" = "energy1",
      "40" = "energy2",
      paste0("energy_", this_ce)
    )
    readr::write_lines(energy.type, file = file, append = TRUE)
    peaks <- Spectra::peaksData(sp.this)
    rows <- lapply(seq_along(peaks), function(i) {
      m <- as.matrix(peaks[[i]])
      if (!nrow(m)) return(character(0))
      o <- order(m[, 1])
      paste0(m[o, 1], " ", m[o, 2])
    })
    readr::write_lines(unlist(rows, use.names = FALSE), file = file, append = TRUE)
  }
  invisible(file)
}

#' @describeIn CFM predict and annotate
#' @description
#' Combines \code{\link{CFM_predict}} and \code{\link{CFM_annotate}} in a single workflow.
#' First predicts mass spectra using CFM-ID, then annotates the predicted spectra
#' with fragment assignments. This is useful for obtaining fully annotated spectral data
#' for a given molecule.
#'
#' Optional disk cache (\code{check_cache}): stores/loads
#' \code{<id>_positive_CFM_data.rds} / \code{<id>_negative_CFM_data.rds}
#' under \code{cache_dir}.
#'
#' @param smiles_or_inchi SMILES string or InChI string of the molecule
#' @param id Compound identifier. If \code{NULL} or empty, a random id is generated
#'   each call (cache will not be reusable across runs unless \code{id} is supplied).
#' @param ppm_mass_tol Mass tolerance in ppm for annotation matching (default: 5.0)
#' @param abs_mass_tol Absolute mass tolerance in m/z units (default: 0)
#' @param param_adduct Adduct type, e.g., "\[M+H\]+" or "\[M-H\]-" (default: "\[M+H\]+").
#'   Also accepts \code{0}/\code{1} (0 = "\[M-H\]-", 1 = "\[M+H\]+") and \code{"-"}/\code{"+"}.
#' @param output_file Path to save results, or NULL to return results in memory (default: NULL)
#' @param check_cache Logical; if TRUE, reuse/save RDS cache under \code{cache_dir}
#' @param cache_dir Cache directory (used when \code{check_cache} is TRUE).
#'   Default \code{tempdir()}.
#' @param force Logical; if TRUE, ignore existing cache and recompute
#' @param ... Additional arguments passed to underlying functions
#'
#' @export
CFM_annotate_by_predict <- function(
    smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
    id = NULL,
    ppm_mass_tol = 5.0,
    abs_mass_tol = 0,
    param_adduct = "[M+H]+",
    output_file = NULL,
    check_cache = FALSE,
    cache_dir = tempdir(),
    force = FALSE,
    ...){

  ### normalize param_adduct (0/1 -> [M-H]- / [M+H]+)
  {
    param_adduct <- .normalize_cfm_param_adduct(param_adduct)
  }

  ### resolve id (random if not supplied)
  {
    if (is.null(id) || !nzchar(as.character(id)[[1]]) || is.na(id)) {
      id <- paste0(
        "CFM_",
        format(Sys.time(), "%Y%m%d%H%M%OS3"),
        "_",
        as.integer(runif(1, 1e8, 1e9 - 1))
      )
      message("[MSCC] CFM id not supplied; using random id: ", id)
    } else {
      id <- as.character(id)[[1]]
    }
  }

  ### check_cache: resolve cache path
  {
    .cfm_cache_polarity_label <- function(adduct) {
      pol <- tryCatch(get_polarity_from_adduct(adduct), error = function(e) NA_integer_)
      if (identical(pol, 0L) || identical(pol, 0)) return("negative")
      if (identical(pol, 1L) || identical(pol, 1)) return("positive")
      # Fallback from common adduct strings
      ad <- as.character(adduct)[[1]]
      if (grepl("\\[M-H\\]-", ad, fixed = TRUE) || identical(ad, "-") || identical(ad, "0")) {
        return("negative")
      }
      "positive"
    }
    .cfm_cache_path <- function(cache_dir, id, adduct) {
      # Polarity label keeps [M+H]+ / [M-H]- distinct without MpHp-style tokens.
      file.path(cache_dir, paste0(id, "_", .cfm_cache_polarity_label(adduct), "_CFM_data.rds"))
    }
    .cfm_cache_candidates <- function(cache_dir, id, adduct) {
      # Prefer new positive/negative names; also accept legacy filenames.
      ad_legacy <- as.character(adduct)
      ad_legacy <- gsub("\\[|\\]", "", ad_legacy)
      ad_legacy <- gsub("\\+", "p", ad_legacy)
      ad_legacy <- gsub("-", "m", ad_legacy)
      ad_legacy <- gsub("[^A-Za-z0-9_]+", "_", ad_legacy)
      unique(c(
        .cfm_cache_path(cache_dir, id, adduct),
        file.path(cache_dir, paste0(id, "_", ad_legacy, "_CFM_data.rds")),
        file.path(cache_dir, paste0(id, "_", adduct, "_CFM_data.rds"))
      ))
    }
  }

  ### check_cache: try load RDS
  {
    cache_file <- .cfm_cache_path(cache_dir, id, param_adduct)
    if (isTRUE(check_cache) && !isTRUE(force)) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      }
      hits <- .cfm_cache_candidates(cache_dir, id, param_adduct)
      hits <- hits[file.exists(hits)]
      if (length(hits)) {
        cached <- tryCatch(readRDS(hits[[1]]), error = function(e) NULL)
        if (inherits(cached, "CFM_data")) {
          message("[MSCC] CFM cache hit: ", normalizePath(hits[[1]], winslash = "/", mustWork = FALSE))
          return(cached)
        }
      }
    }
  }

  ### predict + annotate
  {
    message("[MSCC] CFM_annotate_by_predict: ", id, " ", param_adduct)
    cfm.pred <- CFM_predict(smiles_or_inchi,
                              param_adduct = param_adduct)
    cfm.pred.sp <- get_CFM_data_Spectra(cfm.pred)
    cfm.anno <- CFM_annotate(smiles_or_inchi = smiles_or_inchi, id = id,
                 spectrum_file = cfm.pred.sp,
                 ppm_mass_tol = ppm_mass_tol,
                 abs_mass_tol = abs_mass_tol,
                 param_adduct = param_adduct)
  }

  ### check_cache: save RDS
  {
    if (isTRUE(check_cache) && inherits(cfm.anno, "CFM_data")) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      }
      ok <- tryCatch({
        saveRDS(cfm.anno, cache_file)
        TRUE
      }, error = function(e) {
        warning("CFM cache save failed for ", id, " ", param_adduct, ": ",
                conditionMessage(e), call. = FALSE)
        FALSE
      })
      if (isTRUE(ok)) {
        message("[MSCC] CFM cached: ", normalizePath(cache_file, winslash = "/", mustWork = FALSE))
      }
    }
  }

  return(cfm.anno)
}


#' @describeIn CFM fraggen
CFM_fraggen <- function(smiles_or_inchi = "[H]C1(O)O[C@]([H])(CO)[C@@]([H])(O)[C@]([H])(O)[C@@]1([H])O",
                        max_depth = 2,
                        param_adduct = "[M+H]+",
                        output_file_or_dir = NULL){

  param_adduct <- .normalize_cfm_param_adduct(param_adduct)

  out.file <- ifelse(is.null(output_file_or_dir),
                     tempfile(),
                     normalizePath(output_file_or_dir,mustWork = F))
  CFM(result_dir = dirname(out.file),
      cfmid_cmd = paste0(
        "fraggraph-gen ",
        smiles_or_inchi," ",
        max_depth," ",
        switch(param_adduct,
               "[M+H]+"="+",
               "[M-H]-"="-"),
        " fullgraph ",
        basename(out.file), " "

      ))->cmd

  #message("Make Command to docker: ")
  #message(crayon::red(cmd)%>%crayon::reset() )
  shell(cmd)
  if (is.null(output_file_or_dir)) {
    return(invisible(read_CFM_fraggen_result(out.file)))
  }else{
    #message("Result export to: \n",crayon::red(output_file_or_dir))
    return(invisible(out.file))
  }


}

read_CFM_annotate_result <- function(result_path = "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_annotate_result.txt",
                                     polarity = NA_real_)
  {

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep)
                  )%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = dplyr::case_when(session.type!=0~NA,
                                              T~session.energy)
                  )

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak_assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data),
          d1 =dplyr::case_when(
            assigned~str_extract(line.data,".*(?= \\()"),
            T~line.data
          ),
          fragment_score = str_extract(line.data,"(?<= \\()[^\\)]*"))%>%
        dplyr::mutate(
          mz = str_extract(d1,"^([^\\s]*\\s){1}[^\\s]*"),
          intensity = str_extract(mz,"(?<=[^\\s]{1,10}\\s).*"),
          mz = str_extract(mz,".*(?=[^\\s]{1,10}\\s)"),
          fragment_id = str_extract(d1,"(?<=^[^\\s]{1,10}\\s[^\\s]{1,10}\\s)[^\\s].*"),
          fragment_id = strsplit(x = fragment_id,split = "\\s"),
          fragment_score = strsplit(x = fragment_score,split = "\\s")
        )

      if (nrow(energy.data) == 0) {
        next
      }
      for (j in 1:nrow(energy.data)) {
        this.peak.assign <- data.frame(energy = this.energy,
                        mz = energy.data$mz[j],
                        intensity = energy.data$intensity[j],
                        fragment_id = energy.data$fragment_id[[j]],
                        fragment_score = energy.data$fragment_score[[j]])

        peak_assignment <- bind_rows(peak_assignment,this.peak.assign)

      }


    }
    peak_assignment <- peak_assignment%>%
      dplyr::mutate(mz = as.numeric(mz),
                    intensity = as.numeric(intensity),
                    fragment_score = as.numeric(fragment_score))
    }

  ### session 2
  {
    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
                    )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)
  }

  ### session 3
  {

    session.data.3 <- cfm.df%>%
      dplyr::filter(session.type==2,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)

  }


  ### session 4 fragment transition
  {

    session.data.4 <- cfm.df%>%
      dplyr::filter(session.type==3,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,smiles)%>%
      dplyr::distinct(from,to,smiles)

  }

  ### re id and add polarity suffix
  {
    polarity_suffix <- get_polarity_suffix(polarity)
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.3)), polarity_suffix)%>%
      `names<-`(session.data.3$fragment_id)
    peak_assignment <-peak_assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.3 <- session.data.3%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)
    session.data.4 <- session.data.4%>%
      dplyr::mutate(from = unname(fragment.id.new[from]),
                    to = unname(fragment.id.new[to]))

  }

  ### add polarity column
  {
    peak_assignment$polarity <- polarity
    session.data.3$polarity <- polarity
    session.data.4$polarity <- polarity
  }

  cfm_data <- list(
    peak_assignment = peak_assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.3,
    fragment_transition = session.data.4
  )

  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@peak_assignment <- peak_assignment
    cfm_data@fragment_define <- session.data.3
    cfm_data@fragment_transition <- session.data.4

  }

  return(cfm_data)
}

#' Read CFM-ID Prediction Results
#' @title Read CFM-ID Prediction Results
#' @description Reads and parses the output from CFM-ID prediction results into a structured CFM_data object.
#' This function processes the text output file from CFM-ID and extracts peak assignments,
#' fragment definitions, and other spectral data.
#' @describeIn CFM read prediction results
#' @param result_path Path to the CFM-ID prediction result file (default: "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_predict_result.txt")
#' @param polarity Numeric polarity value (0 for negative, 1 for positive, default: NA)
#'
#' @return A CFM_data object containing parsed peak assignments and fragment definitions
#' @export
#'

read_CFM_predict_result <- function(result_path =  "c:/Users/91879/OneDrive/Code/Docker/cfm/data/cfm_predict_result.txt",
                                    polarity = NA_real_){

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep))
  cfm.df <- cfm.df[!grepl(pattern = "^#",  x = cfm.df$line.data),]
    #dplyr::mutate(
    #  session.type = dplyr::case_when(grepl(pattern = "^#",
    #                                             x = line.data)~-1,
    #                                       T~session.type)
    #)%>%
  cfm.df <- cfm.df%>%
    dplyr::mutate(session.energy.sep = grepl(pattern = "energy",x = line.data ),
                  session.energy = cumsum(session.energy.sep),
                  session.energy = dplyr::case_when(session.type!=0~NA,
                                             T~session.energy)
    )
  ### header
  {


    }

  ### session 1
  {
    energy <- unique(cfm.df$session.energy)%>%na.omit()
    peak_assignment <- data.frame()
    for (i in seq_along(energy)  ) {
      this.energy <-cfm.df%>%
        dplyr::filter(session.energy==energy[i],session.energy.sep)%>%
        dplyr::pull(line.data)
      energy.data <- cfm.df%>%
        dplyr::filter(session.energy==energy[i],!session.energy.sep)%>%
        dplyr::mutate(
          assigned = grepl(pattern = "\\(",x = line.data))%>%
        dplyr::mutate(
          d1 = case_when(
            assigned~str_extract(line.data,".*(?= \\()"),
            T~line.data
          ),
          fragment_score = str_extract(line.data,"(?<= \\()[^\\)]*"))%>%
        dplyr::mutate(
          mz = str_extract(d1,"^([^\\s]*\\s){1}[^\\s]*"),
          intensity = str_extract(mz,"(?<=[^\\s]{1,10}\\s).*"),
          mz = str_extract(mz,".*(?=[^\\s]{1,10}\\s)"),
          fragment_id = str_extract(d1,"(?<=^[^\\s]{1,10}\\s[^\\s]{1,10}\\s)[^\\s].*"),
          fragment_id = strsplit(x = fragment_id,split = "\\s"),
          fragment_score = strsplit(x = fragment_score,split = "\\s")
        )

      if (nrow(energy.data) == 0) {
        next
      }
      for (j in 1:nrow(energy.data)) {
        this.peak.assign <- data.frame(energy = this.energy,
                                       mz = energy.data$mz[j],
                                       intensity = energy.data$intensity[j],
                                       fragment_id = energy.data$fragment_id[[j]],
                                       fragment_score = energy.data$fragment_score[[j]])

        peak_assignment <- bind_rows(peak_assignment,this.peak.assign)

      }


    }
    peak_assignment <- peak_assignment%>%
      dplyr::mutate(mz = as.numeric(mz),
                    intensity = as.numeric(intensity),
                    fragment_score = as.numeric(fragment_score))
  }

  ### session 2
  {

    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles)

  }

  ### re id and add polarity suffix
  {
    polarity_suffix <- get_polarity_suffix(polarity)
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.2)), polarity_suffix)%>%
      `names<-`(session.data.2$fragment_id)
    peak_assignment <-peak_assignment%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
    session.data.2 <- session.data.2%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)

  }

  ### add polarity column
  {
    peak_assignment$polarity <- polarity
    session.data.2$polarity <- polarity
  }

  cfm_data <- list(
    peak_assignment = peak_assignment,
    #fragment_define1 = session.data.2,
    fragment_define = session.data.2
  )


  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@peak_assignment <- peak_assignment
    cfm_data@fragment_define <- session.data.2


  }
  return(invisible(cfm_data))

}

read_CFM_fraggen_result <- function(result_path){

  cfm_data <- readr::read_lines(result_path)
  cfm.df <- data.frame(line.no = 1:length(cfm_data),
                       line.data = cfm_data)%>%
    dplyr::mutate(session.sep = line.data == "",
                  session.type = cumsum(session.sep),
                  session.type = dplyr::case_when(grepl(pattern = "^#",
                                                 x = line.data)~-1,
                                           T~session.type)
    )

  ### session 1 fragment def
  {

    session.data.1 <- cfm.df%>%
      dplyr::filter(session.type==0,!session.sep)%>%
      dplyr::mutate(intermediate = grepl(pattern =" Intermediate Fragment",x = line.data),
                    line.data= gsub(pattern =" Intermediate Fragment",x = line.data,replacement= ""))%>%
      dplyr::mutate(fragment_id  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    fragment_mz = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    fragment_mz = as.numeric(fragment_mz),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(fragment_id))%>%
      dplyr::select(fragment_id,fragment_mz,smiles,intermediate)

    }


  ### session 2 fragment transition
  {

    session.data.2 <- cfm.df%>%
      dplyr::filter(session.type==1,!session.sep)%>%
      dplyr::mutate(from  = str_extract(line.data,"^[^\\s]*(?=\\s)"),
                    to = str_extract(line.data,"(?<=\\s)[^\\s]*(?=\\s)"),
                    neutral_loss = str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s]*(?=\\s)"),
                    smiles =  str_extract(line.data,"(?<=^[^\\s]{1,50}\\s[^\\s]{1,50}\\s[^\\s]{1,50}\\s)[^\\s].*")
      )%>%
      dplyr::filter(!is.na(from))%>%
      dplyr::select(from,to,neutral_loss,smiles)

  }

  ### re id
  {
    fragment.id.new <- paste0("Fragment",num2str(1:nrow(session.data.1)))%>%
      `names<-`(session.data.1$fragment_id)
    session.data.1 <- session.data.1%>%
      dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id] ))%>%
      `rownames<-`(.$fragment_id)
    session.data.2 <- session.data.2%>%
      dplyr::mutate(from = unname(fragment.id.new[from]),
                    to = unname(fragment.id.new[to]))

  }


  ### create cfm_data
  {
    cfm_data <- new("CFM_data")
    cfm_data@fragment_define <- session.data.1
    cfm_data@fragment_transition <- session.data.2

  }


  return(invisible(cfm_data))

}



plot_CFM_annotated_Spectra <- function(cfmd){

  peak.assign <- cfmd@peak_assignment
  fragment.define <- cfmd@fragment_define

  get_CFM_data_Spectra(cfmd) %>%
    Spectra::combineSpectra() %>%
    MSdev::plot_Spectra(label.top = 0)

}



#' Annotate Isotopologues in Mass Spectra
#' @title Annotate Isotopologues in Mass Spectra
#' @description Annotates mass spectra with isotopologue information by matching peaks
#' to expected isotope patterns. This function assigns fragment groups and isotope counts
#' to observed peaks based on mass tolerance.
#' @describeIn MSIP annotate isotopologues
#'
#' @param sp A Spectra object containing experimental mass spectra data
#' @param cfmd A CFM_data object containing fragment definitions and peak assignments
#' @param iso_ele Isotope element specification, e.g., "\[13\]C" for carbon-13 (default: "\[13\]C")
#' @param iso_count Maximum number of isotope incorporations to consider (default: 0)
#' @param ppm Mass tolerance in parts per million for isotope matching (default: 20)
#'
#' @return A data frame containing annotated spectrum data with fragment groups, isotope counts,
#' intensity ratios, and summed intensities for each fragment group
#' @export
#'
CFM_annotate_isotopologues <- function(sp,
                                 cfmd,
                                 iso_ele = "[13]C",
                                 iso_count = 0 ,
                                 ppm = 20){

  # Check if fragment_group exists in peak_assignment
  # For MSIPAtomMap, fragment_group should already be set by MSIPAtomMap_get_FG_map
  if (!"fragment_group" %in% colnames(cfmd@peak_assignment)) {
    # If fragment_group is not in peak_assignment, try to match from fragment_define
    if ("fragment_group" %in% colnames(cfmd@fragment_define)) {
      cfmd@peak_assignment$fragment_group <-
        cfmd@fragment_define$fragment_group[match(
          cfmd@peak_assignment$fragment_id,
          cfmd@fragment_define$fragment_id)]
    } else {
      warning("fragment_group not found. Run MSIPAtomMap_get_FG_map() first.")
    }
  }

  cfm.peaks.data <- cfmd@peak_assignment%>%
    dplyr::filter(!is.na(fragment_id))%>%
    dplyr::mutate(collisionEnergy = case_when(energy == "energy0"~10,
                                              energy == "energy1"~20,
                                              energy == "energy2"~40,
    ))
  if (!nrow(cfm.peaks.data))
    cfm.peaks.data <-  cfmd@peak_assignment%>%
    dplyr::mutate(mz = 0)

  diff.formula <- paste0(iso_ele,MSCC::get_ele_uniso(iso_ele),"-1")
  iso.mz.diff <- (0:iso_count)*MSCC::chemform_mz(diff.formula)
  mz.labeled.m <- matrixSub(cfm.peaks.data$mz,-iso.mz.diff)%>%
    `colnames<-`(paste0("M",0:iso_count))%>%
    as.data.frame()%>%
    dplyr::mutate(fragment_group=cfm.peaks.data$fragment_group)%>%
    dplyr::distinct(M0,.keep_all = T)%>%
    tidyr::pivot_longer(paste0("M",0:iso_count),
                        names_to = "iso_count",
                        values_to = "mz")



  sp.data <- get_Spectra_data(sp,var = "collisionEnergy")%>%
    dplyr::mutate(fragment_group = NA,iso_count=NA)
  if (nrow(sp.data)) {


  sp.data <- sp.data%>%
    dplyr::mutate(idx = match_mz(mz1=mz,
                                mz2 = mz.labeled.m$mz,
                                mz.ppm = ppm),
                  fragment_group = mz.labeled.m$fragment_group[idx],
                  iso_count = mz.labeled.m$iso_count[idx])%>%
    dplyr::select(-idx)%>%
    dplyr::group_by(fragment_group,sp.id)%>%
    dplyr::mutate(
      ratio = case_when(
      !is.na(fragment_group)~intensity/sum(intensity),
      T~NA),
      int_sum =  case_when(
        !is.na(fragment_group)~sum(intensity),
        T~NA))%>%
    dplyr::ungroup()
  }


  return(sp.data)
}


#' Get Igraph Objects for CFM Fragments
#' @title Get Igraph Objects for CFM Fragments
#' @description Converts fragment structures in a CFM_data object to igraph objects representing
#' molecular graphs. This function processes fragment SMILES strings to SDF format and then
#' creates igraph objects for structural analysis.
heatmap_atom_iso_prob <- function(x){

  ComplexHeatmap::Heatmap(x,
                          na_col  ="#999999",
                          name = "isotope labeled\nprobability",
                          col = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                     c("white","#F7844F","#B20C26")),
                          cluster_columns = F,
                          row_names_side  = "left",
                          rect_gp =  grid::gpar(lwd=2,col = "white"),
                          cluster_rows = F)

}


#' Weight CFM Spectra Intensities by Fragment Group
#' @title Weight CFM Spectra Intensities by Fragment Group
#' @description Computes intensity-weighted mean isotopologue intensities for each
#' fragment group across spectra, then appends a combined spectrum
#' (\code{sp.id = "combined_sp"}) to the input data.
#'
#' @param sp.data Data frame of annotated spectrum peaks, typically from
#'   \code{\link{CFM_annotate_isotopologues}}, with columns such as
#'   \code{fragment_group}, \code{sp.id}, \code{iso}, \code{intensity}, and \code{mz}.
#' @param iso_count Maximum isotope count (Mn) to include, as an integer.
#'
#' @return A data frame with the original peaks plus weighted combined peaks.
#' @export
CFM_spectra_data_int_weight <- function(sp.data,iso_count){



  ### weighted intensity matrix
  {
    fg.idx <- split(1:nrow(sp.data),sp.data$fragment_group)
    if (!length(fg.idx))  return(sp.data)
    frag.iso.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    frag.int.sum <- c()
    for (i.fg in seq_along(fg.idx)) {
      x.df <- sp.data[fg.idx[[i.fg]],]
      x.int <- x.df%>%
        dplyr::select(-mz,-collisionEnergy)%>%
        tidyr::pivot_wider(names_from ="iso",
                           id_cols = "sp.id",
                           values_from = "intensity",
                           values_fn = sum)%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      to.add <- setdiff(paste0("M",0:iso_count),colnames(x.int))
      x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                            dimnames = list(NULL,to.add)),x.int)
      x.int <- x.int[,paste0("M",0:iso_count),drop =F]
      x.int[is.na(x.int)] <- 0
      if (ncol(x.int)==1) {
        x.ratio <- x.int
        x.ratio[,1] <- 1
      }else{
        x.ratio <- t(apply(x.int,1,function(z) z/sum(z)))
      }
      x.weight <- rowSums(x.int)
      x.int.weighted <- apply(x.int,2,weighted.mean,w = x.weight)
      x.ratio.weighted <- apply(x.ratio,2,weighted.mean,w = x.weight)
      frag.iso.matrix[i.fg,] <- x.int.weighted
      frag.int.sum[i.fg] <- sum(x.int.weighted)
    }
    names(frag.int.sum) <- names(fg.idx)

  }


  ### merge
  {
    sp.mz <- sp.data%>%
      dplyr::filter(!is.na(fragment_group))%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,iso))

    iso.peak.data <- frag.iso.matrix%>%
      as.data.frame()%>%
      rownames_to_column("fragment_group")%>%
      tidyr::pivot_longer(starts_with("M"),
                          names_to = "iso",
                          values_to = "intensity")%>%
      dplyr::mutate(sp.id= "combined_sp",
                    fg.iso = paste0(fragment_group,iso),
                    mz = sp.mz$mz[match(fg.iso,sp.mz$fg.iso)])%>%
      dplyr::filter(!is.na(mz))
    sp.data.weighted <- bind_rows(
      sp.data,
      iso.peak.data
    )%>%
      dplyr::select(-fg.iso)


  }
  return(sp.data.weighted)
}

#' Merge CFM Spectra Isotopologue Ratios by Fragment Group
#' @title Merge CFM Spectra Isotopologue Ratios by Fragment Group
#' @description Merges per-spectrum isotopologue ratios within each fragment group
#' using intensity-weighted means, and appends combined peaks
#' (\code{sp.id = "combined_sp"}, \code{merged = TRUE}) with summary metrics
#' (\code{int_sum}, \code{peaks_count}, \code{icc}, \code{cos}).
#'
#' @param sp.data Data frame of annotated spectrum peaks with fragment-group
#'   ratios, typically from \code{\link{CFM_annotate_isotopologues}}, including
#'   columns such as \code{fragment_group}, \code{sp.id}, \code{iso_count},
#'   \code{ratio}, \code{int_sum}, and \code{mz}.
#' @param iso_count Maximum isotope count (Mn) to include, as an integer.
#'
#' @return A data frame with original peaks plus merged combined peaks and metrics.
#' @export
CFM_spectra_data_merge <- function(sp.data,iso_count){


  sp.data$merged <- F
  ### weighted intensity matrix
  {
    fg.idx <- split(1:nrow(sp.data),sp.data$fragment_group)
    if (!length(fg.idx))  return(sp.data)
    frag.ratio.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    frag.df <- data.frame(
      fragment_group = names(fg.idx),
      int_sum = NA,
      peaks_count = NA,
      icc = NA,
      cos = NA
    )
    for (i.fg in seq_along(fg.idx)) {

      x.df <- sp.data[fg.idx[[i.fg]],]
      x.ratio <- x.df%>%
        tidyr::pivot_wider(names_from ="iso_count",
                           id_cols = "sp.id",
                           values_from = "ratio",
                           values_fn = mean
                           )%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      x.ratio <- get_matrix_value_fill_with_NA(
        x.ratio,
        rownames_vec = row.names(x.ratio),
        colnames_vec = paste0("M",0:iso_count),drop = F)
      x.ratio[is.na(x.ratio)] <- 0
      if (ncol(x.ratio)==1) {
        x.ratio[,1] <- 1
      }
      x.int.sum <- x.df%>%
        dplyr::distinct(sp.id,.keep_all = T)%>%
        dplyr::pull(int_sum,name = sp.id)
      x.int.sum <- x.int.sum[rownames(x.ratio)]
      x.ratio.weighted <- apply(x.ratio,2,weighted.mean,w = x.int.sum)
      x.int.sum.weighted <- weighted.mean(x.int.sum,x.int.sum)
      x.icc <- irr::icc(t(x.ratio), model = "twoway",
          type = "consistency", unit = "single")$value
      x.cos <- lsa::cosine(x.ratio.weighted,t(x.ratio))
      x.cos.weight <- weighted.mean(x.cos,w = log10(x.int.sum))
      frag.ratio.matrix[i.fg,] <- x.ratio.weighted
      frag.df$int_sum[i.fg] <- x.int.sum.weighted
      frag.df$icc[i.fg] <- x.icc
      frag.df$cos[i.fg] <- x.cos.weight
      frag.df$peaks_count[i.fg] <- nrow(x.ratio)

      #Heatmap(x.ratio,cluster_columns = F,cluster_rows = F)

    }

  }


  ### merge
  {
    sp.mz <- sp.data%>%
      dplyr::filter(!is.na(fragment_group))%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,iso_count))

    iso.peak.data <- frag.ratio.matrix%>%
      as.data.frame()%>%
      rownames_to_column("fragment_group")%>%
      tidyr::pivot_longer(starts_with("M"),
                          names_to = "iso_count",
                          values_to = "ratio")%>%
      dplyr::mutate(sp.id= "combined_sp",
                    fg.iso = paste0(fragment_group,iso_count),
                    mz = sp.mz$mz[match(fg.iso,sp.mz$fg.iso)],
                    frag.df[match(fragment_group,frag.df$fragment_group),],
                    intensity = int_sum*ratio ,
                    merged = T )%>%
      dplyr::filter(!is.na(mz))

    sp.data$merged <- F
    sp.data.weighted <- bind_rows(
        sp.data,
      iso.peak.data )%>%
      dplyr::select(-fg.iso)


  }

  return(sp.data.weighted)
}

#' Remove Natural Isotope Contribution from Combined CFM Spectra
#' @title Remove Natural Isotope Contribution from Combined CFM Spectra
#' @description Subtracts estimated natural isotopologue intensity from combined
#' spectrum peaks (\code{sp.id = "combined_sp"}) using an isoform map and a
#' natural-abundance scaling factor. Peaks with non-positive intensity after
#' subtraction are set to zero; fragment groups that become all-zero are dropped
#' from the combined spectrum.
#'
#' @param sp.data Data frame of spectrum peaks that may include a combined
#'   spectrum from \code{\link{CFM_spectra_data_int_weight}} or
#'   \code{\link{CFM_spectra_data_merge}}.
#' @param natural.ratio Numeric scaling factor for the natural isotope contribution.
#' @param if.map Object with an \code{isoform.map} slot used to derive the natural
#'   isotopologue distribution.
#'
#' @return A data frame with natural contribution removed from combined peaks.
#' @export
CFM_spectra_data_remove_natural <-function(sp.data,
                                           natural.ratio,
                                           if.map){

  ###
  {
    if (sum(sp.data$sp.id=="combined_sp")==0) return(sp.data)
  }


  ### calculate intensity from natrual
  {
    sp.data.merged <- sp.data%>%
      dplyr::filter(sp.id== "combined_sp")%>%
      dplyr::mutate(fg.iso = paste0(fragment_group,"_",iso))
    int_sum <- sp.data.merged%>%
      dplyr::group_by(fragment_group)%>%
      dplyr::summarise(int_sum = sum(intensity))%>%
      #dplyr::mutate(int_sum = sum(intensity),
      #              fgi = paste0(fragment_group,"_",iso))%>%
      #column_to_rownames("fragment_group")%>%
      dplyr::pull(int_sum,name  = fragment_group)
    max.iso <- max(str_extract_num(sp.data.merged$iso))
    suffix <- paste0("_M",rep(0:max.iso,length(int_sum)))
    int_sum <- rep(int_sum,each = max.iso+1)
    names(int_sum) <- paste0(names(int_sum) , suffix)
    natural.distribution <- apply(if.map@isoform.map,1,
                                sum)/ncol(if.map@isoform.map)
    natural.int <- int_sum[names(natural.distribution)]*natural.distribution*natural.ratio

  }

  {
    sp.data.removed <-  sp.data.merged %>%
      dplyr::mutate(intensity = intensity - natural.int[fg.iso],
                    intensity = case_when(intensity <0 ~0,
                                          is.na(intensity)~0,
                                          T~intensity))%>%
      dplyr::group_by(fragment_group)%>%
      dplyr::filter(!all(intensity==0))%>%
      dplyr::select(-fg.iso)
    sp.data <- sp.data%>%
      dplyr::filter(!sp.id== "combined_sp")%>%
      rbind(sp.data.removed)
  }

  return(sp.data)


}


get_CFM_data_MSIPFragmentMap<- function(msipAtomMap){

  cfmd.sp <- get_CFM_data_Spectra(msipAtomMap)
  cfmd.msip.core <- get_MSIPCoreData(cfmd.sp,msipAtomMap,0)
  cfmd.fg.map <- cfmd.msip.core@MSIPFragmentMap

  cfmd.fg.map
}



#' Get CFM Data from SMILES
#' @title Get CFM Data from SMILES
#' @description Creates CFM_data object from SMILES by running CFM prediction and annotation.
#'
#' @param smiles SMILES string of the molecule
#' @param compound_id Identifier for the compound. If \code{NULL} or empty, a
#'   random id is generated each call (see \code{\link{CFM_annotate_by_predict}}).
#' @param ppm Mass tolerance in ppm (default: 5)
#' @param adduct Adduct type (default: "\[M+H\]+"). Also accepts 0/1 (0 = "\[M-H\]-", 1 = "\[M+H\]+")
#'   and '-'/'+' as shorthand.
#' @param check_cache Logical; forward to \code{\link{CFM_annotate_by_predict}}
#' @param cache_dir Cache directory; forward to \code{\link{CFM_annotate_by_predict}}
#'   (default \code{tempdir()}).
#' @param force Logical; ignore cache if TRUE
#' @param ... Additional arguments
#'
#' @return A CFM_data object containing fragment data
#' @export
get_CFM_data_from_smiles <- function(smiles = "NCC(O)=O",
                                     compound_id = NULL,
                                     ppm = 5,
                                     adduct = "[M+H]+",
                                     check_cache = FALSE,
                                     cache_dir = tempdir(),
                                     force = FALSE,
                                     ...){
  adduct <- .normalize_cfm_param_adduct(adduct)

  cfm_data <- CFM_annotate_by_predict(
    smiles_or_inchi = smiles,
    id = compound_id,
    ppm_mass_tol = ppm,
    abs_mass_tol = 0.005,
    param_adduct = adduct,
    check_cache = check_cache,
    cache_dir = cache_dir,
    force = force
  )
  return(cfm_data)
}



shiny_vis_cfmd_FG_map <- function(msipAtomMap){


  .ui <- function(){
    fluidPage(
      column(width = 6,
             shiny::plotOutput(outputId = "heatmap_atom_map",height  = "800px")),
      column(width = 6,
             selectInput(inputId = "fg_id",choices = sort(unique(msipAtomMap@fragment_define$fragment_group)),
                         label = "Fragment group",selected = msipAtomMap@fragment_define$fragment_group[1]),
             selectInput(inputId = "fragment",label = "fragment",
                         choices = sort(unique(msipAtomMap@fragment_define$fragment_id)),
                         selected = msipAtomMap@fragment_define$fragment_id[1] ),
             visNetwork::visNetworkOutput(outputId = "atom_map",height  = "800px"))
    )

  }

  .server <- function(){
    function(input, output, session) {

      output$heatmap_atom_map <- renderPlot({

        message_with_time("heatmap_atom_map")
        get_MSIPAtomMap_MSIPFragmentMap(msipAtomMap)%>%
          heatmap_MSIPFragmentMap()
      })

      output$atom_map <- visNetwork::renderVisNetwork({
        message_with_time("atom_map")
        vis_MSIPAtomMap_fragment_atom_map(msipAtomMap ,input$fragment,show_id = F)
      })

      observeEvent(input$fg_id,{

        message_with_time("fg_id")

        x <- msipAtomMap@fragment_define%>%
          dplyr::filter(fragment_group == input$fg_id)
        updateSelectInput(inputId = "fragment",choices = x$fragment_id)
      })

    }
  }

  ### Start Shiny APP
  {
    shinyApp(ui = .ui(),
             server = .server(),
             options = list(host = "0.0.0.0",
                            #port = 6548,
                            launch.browser = T))
  }


}


shiny_vis_cfmd_trans <- function(msipAtomMap){


  .ui <- function(){
    fluidPage(
      column(width = 6,
             visNetwork::visNetworkOutput(outputId = "vis_trans",height = "800px"),
             style = "border: 1px solid #aaa; padding: 6px;"),
      column(width = 6,
             plotlyOutput(outputId = "cfm_sp"),
             verbatimTextOutput("frag_info"),
             visNetwork::visNetworkOutput(outputId = "atom_map" ))
    )

  }

  .server <- function(){
    function(input, output, session) {

      output$vis_trans <-  visNetwork::renderVisNetwork({

        message_with_time("get_MSIPAtomMap_trans_igraph")
        vis_igraph(get_MSIPAtomMap_trans_igraph(msipAtomMap)) %>%
          visOptions(nodesIdSelection =
                       list(enabled  = T,
                            selected  = msipAtomMap@fragment_define$fragment_id[1]))


      })


      output$frag_info <- renderText({

        message_with_time("frag_info")
        x <- input$vis_trans_selected
        if(is.null(x)) return(NULL)
        if(x==""){      return(NULL)  }
        y <- msipAtomMap@fragment_define
        paste0("Fragment: ",x,"\n",
               "Formula: ",y[x,"formula"],"\n",
               "mz: ",y[x,"fragment_mz"],"\n"  )

      })
      output$atom_map <-  visNetwork::renderVisNetwork({

        message_with_time("get_MSIPAtomMap_trans_igraph")
        x <- input$vis_trans_selected
        if(is.null(x)) return(NULL)
        if(x==""){      return(NULL)  }

        vis_MSIPAtomMap_fragment_atom_map(msipAtomMap,input$vis_trans_selected,show_id = F)


      })


      output$cfm_sp <- renderPlotly({

        MSdev::plotly_Spectra(get_CFM_data_Spectra(msipAtomMap))

      })


    }
  }

  ### Start Shiny APP
  {
    shinyApp(ui = .ui(),
             server = .server(),
             options = list(host = "0.0.0.0",
                            #port = 6548,
                            launch.browser = T))
  }


}




get_cfmd_FG_map_check <- function(msipAtomMap){

  fg.map <- get_MSIPAtomMap_MSIPFragmentMap(msipAtomMap)
  fg.fragment.count <- table(msipAtomMap@fragment_define$fragment_group)
  fg.certainty <- get_MSIPFragmentMap_certainty(fg.map)
  fg.df <- data.frame(
    fg = names(fg.certainty),
    certainty =fg.certainty,
    fragment.count = as.numeric(fg.fragment.count[names(fg.certainty)])
  )

  return(fg.df)

}


get_cfmd_map_evaluation <- function(cfmd){



}

