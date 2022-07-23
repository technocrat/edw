receiver <- data.table(
# lq derives from ss, from JDAY.x set by process_pixels_3g_v2_01[1].R
# better to recreate here, but will use object JDAY.x.rds from Edw. M
# JDAY.x is a vector of length 780 integers, applicable to all records
# pix is the record id added corresponding to original row order
  pix = as.integer(),
  # variables derived from the smoothed splines
  lq = as.double(),
  uq = as.double(),
  mean.evi = as.double(),
  sd.evi = as.double(),
  sum.evi = as.double(),
  amplitude = as.double(),
  # variables derived from the smoother splines and derivatives
  peak.day = as.double(),
  trough.day = as.double(),
  devi.ss = as.double(),
  devi.sss = as.double(),
  # a constant from cons.R
  L.JDAY = as.integer(),
  # variables derived from from scaled data
  NA.length = as.integer(),
  ICE.length = as.integer()
) 

# remaining variables are either undefined
# or supposed to be lists, which will be
# looked at later
