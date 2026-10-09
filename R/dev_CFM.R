plot_CFM_annotated_Spectra <- function(cfmd) {
  MSAtomMap::plotly_CFM_spectra(cfmd)
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
