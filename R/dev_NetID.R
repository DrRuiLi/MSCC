#' @title  Expand compound to adduct
#' @description Using formula, exact.mass and ion_mode, expand adduct by enviPat::adduct
#'
#' @param chem_formula a chemical formula, such as "C4H4Cl1N3O1".
#' some template see "MSCC::chem_formula_template"
#' @param ion_mode
#'
#' @return
#' @export
#'
#' @examples
expand_adduct_from_formula <- function(chem_formula, ion_mode = "positive") {


  #ion_mode <- "positive"
  #chem_formula <- chem_formula_template[1]
  adduct.rule <- adduct.table %>%
    dplyr::filter(Ion_mode == ion_mode)


    x.formula <- chem_formula
    x.mass <- enviPat::check_chemform(isotopes,x.formula)$monoisotopic_mass
    x.isopat <- isotopes_pattern_enviPat(x.formula)%>%
      mutate(form = paste0(isotope_element , "M"))

    x_adduct <-data.frame()
    for (i in 1:nrow(x.isopat)) {

      m_formula <- x.isopat$formula[i]
      m_mass <- x.isopat$m.z[i]
      m_form <- x.isopat$form[i]
      m_adduct <- adduct.rule %>%
        rowwise() %>%
        mutate(
          formula =  enviPat::multiform(formula_in = m_formula , fact = Multi),
          formula = formula_calculate_lc8(formula , Formula_diff),
          adduct = Name,
          adduct = sub(pattern = "M",replacement = m_form , x = adduct),
          charge = Charge,
          multi = Multi,
          ion_mode = Ion_mode,
          form = "adduct",
          exact.mz = m_mass * Multi / abs(Charge) + Mass
        ) %>%
        select("formula"  ,
               "adduct",
               "charge",
               "multi",
               "ion_mode",
               "form",
               "exact.mz")
      m_adduct
      x_adduct <- rbind(x_adduct , m_adduct)

    }



    x_adduct


  return(x_adduct)

}
