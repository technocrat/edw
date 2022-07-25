receiver <- data.table(
# pix is the record id added corresponding to original row order
  pix = as.integer(),
  # variables derived from the smoothed splines
  lq = as.double(),
  uq = as.double(),
  mean.evi = as.double(),
  sd.evi = as.double(),
  sum.evi = as.double(),
  amplitude = as.double(),
  devi.ss = as.double(),
  devi.sss = as.double(),
  # a constant from cons.R
  L.JDAY = as.integer(),
  # variables derived from from scaled data
  NA.length = as.integer(),
  ICE.length = as.integer()
) 

