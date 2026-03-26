


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


