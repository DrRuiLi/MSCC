benchmark <- function(){


  time.stat <- data.frame(number = 1:7,
                          time = NA)
  for (i in 1:7) {
    n <- 10^(i)
    chemform <- MSCC::demo_chemform[sample(1:107,n,replace = T)]
    t <- system.time(
      #chemform_multi(chemform,multi = sample(1,n,replace = T))
     # chemform_calc(chemform,chemform,calc = "-")
      #chemform_parse(chemform,return = "matrix")
      #chemform_add_num(chemform)
      #chemform_formate(chemform,return = "ddd")
      #chemform_mz(chemform)
      chemform_adduct(chemform)
      )
    message(n)
    print(t)
    time.stat$number[i] <- n
    time.stat$time[i] <- t[3]
  }



}
