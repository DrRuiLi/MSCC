


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


  chemform_decompose_mass(203.0526)
  chemform_decompose_mz(115.054399, charge = 1)
  chemform_decompose_mz(117.07004, charge = 1)

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
      elements = c("C", "H", "N", "O", "P")
    )$formula

    stopifnot(setequal(r, mm))
  }

  message("Rdisop vs MSCC parity checks passed (formula sets).")


}

# Wed Jul 29 19:03:00 2026 — seven golden rules filter
{
  f <- c("C6H12O6", "C26H28N17OP3S8", "CH6N2", "C8HN5")
  all_df <- MSCC::chemform_check_seven_golden_rules(f, return = "all")
  stopifnot(identical(
    MSCC::chemform_check_seven_golden_rules(f, return = "valid"),
    all_df$valid
  ))
  stopifnot(identical(
    MSCC::chemform_check_seven_golden_rules(f, return = "chemform"),
    f[all_df$valid]
  ))
  stopifnot(isTRUE(all_df$valid[1]))                 # glucose passes
  stopifnot(all_df$SENIOR3[1] == 2)                  # Senior3 for C6H12O6
  stopifnot(!all_df$elementHeuristic[2])             # high NOPS fails #6
  stopifnot(!all_df$HCratio[3], !all_df$HCratio[4]) # H/C extremes fail #4
  message("chemform_check_seven_golden_rules smoke checks passed.")

  # check_rule on decompose: filtered set is subset of unfiltered
  raw <- MSCC::chemform_decompose_mass(180.0634, ppm = 5, check_rule = FALSE)
  filt <- MSCC::chemform_decompose_mass(180.0634, ppm = 5, check_rule = TRUE)
  stopifnot(nrow(filt) <= nrow(raw))
  stopifnot(all(filt$formula %in% raw$formula))
  if (nrow(filt)) {
    stopifnot(all(MSCC::chemform_check_seven_golden_rules(filt$formula, mass = filt$exactmass, return = "valid")))
  }
  raw_mz <- MSCC::chemform_decompose_mz(181.0707, charge = 1, ppm = 5, check_rule = FALSE)
  filt_mz <- MSCC::chemform_decompose_mz(181.0707, charge = 1, ppm = 5, check_rule = TRUE)
  stopifnot(nrow(filt_mz) <= nrow(raw_mz))
  message("chemform_decompose check_rule smoke checks passed.")
}

{
  target.mz <- 117.07004
  chemform_decompose_mz(target.mz,charge = 1)

}
