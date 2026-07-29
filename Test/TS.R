


time.stat <- c(testn,st)
for (i in 1:6) {

  testn <- 10^i
  chemform <- sample(MSCC::chem_formula_template,testn,replace = T)
  adduct <- sample(MSCC::adduct.table$Adduct,testn,replace = T)
  charge <-sample(-2:2,testn,replace = T)
  #chemform_test.split <- split(adduct_test,
  #                             sample(letters,testn,replace = T))
  st <- system.time(
      chemform_mz(chemform,charge)
  )
  message(st[3])
  time.stat <- rbind(time.stat,c(testn,st))
}



# Wed Jul 29 15:31:06 2026 ------------------------------
{


  library(Rdisop)

  Rdisop::decomposeMass()


  chemform_decompose_mass(masses,charge = 0)

  # Parity checks: formula sets vs Rdisop::decomposeMass(maxisotopes=1)
  els <- initializeElements(c("C", "H", "N", "O", "P"))
  masses <- c(203.0526, 100.0)
  ppm <- 5

  for (m in masses) {
    r <- decomposeMass(
      mass = m,
      ppm = ppm,
      elements = els,
      z = 0,
      maxisotopes = 1
    )$formula

    mm <- MSCC::chemform_decompose_mass(
      mass = m,
      ppm = ppm,
      charge = 0,
      elements = c("C", "H", "N", "O", "P")
    )$formula

    stopifnot(setequal(r, mm))
  }

  message("Rdisop vs MSCC parity checks passed (formula sets).")


}
