


work_flow_mzML <- function(project.dir = "D:/MSCC.test/",
                      compound.formula.file = "d:/MSCC.test/MSCC.compound.xlsx"){



  msdev <- MSdev:::MSdev(rawDataDir = paste0(project.dir,"/rawData"))
  msdev <- MSdev::checkSampleInfo(msdev)
  msdev <- MSdev::msConvert_MSdev(msdev)
  msdev  <- MSdev:::xcmsProcessingMSdev(msdev,
                                        xcms.findpeak.param =xcms::CentWaveParam(ppm = 20,
                                                                                snthresh = 100,
                                                                                peakwidth = c(5, 50),
                                                                                prefilter = c(3, 1000)))



  compound.table <- readxl::read_excel("d:/MSCC.test/MSCC.compound.xlsx")%>%
    dplyr::rowwise()%>%
    dplyr::mutate(Chem_formula = chemform_formate(Chem_formula),
                  Chem_formula_adduct = chemform_adduct(Chem_formula , Adduct)
                  )
  compound.isotopes.network <- BiocParallel::bplapply(X = compound.table$Chem_formula_adduct,
                                                      chemform_isotopes_pattern_enviPat,
                                      BPPARAM = BiocParallel::SerialParam(progressbar = T)
                                                    )

  compound.isotopes.network.match <- match_isotopes_to_features(isotopes.network = compound.isotopes.network,
                                                                xcms.xcms =msdev@xcmsData$positiveMS1,
                                                                ppm.thresh = 20
                                                          )






}
