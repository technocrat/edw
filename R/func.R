# func.R
# functions used in scripts under this directory
# author: Richard Careaga
# Date: 2022-07-20

# takes a row vector of the data and return a scaled vector of the
# same size
# TODO: this needs to be on a row-by-row basis
find_frozen <- function(x){
  # how many possible ICE/SNOW data points
  # remove record id
  x   = x[-c(1)]
  r = x - floor(x/10)*10 +1
  ICE.length = length(which(r==4 | r==6))  
  return(ICE.length)
}

# no NA check is run because approx can be run on rows
# without NAs and makes no changes and is not an expensive
# operation
make_splines <- function(x) {
  # set aside record id
  pix = x[,1]
  x   = x[,-c(1)]
  # get a filtered running mean:
  splines = filter(x,
        # 2022-07-21 00:20 -07:00 
        # convolution won't work with filter of length 3
        filter = 1 / 3, # rep(1/3,3),
        method = "convolution",
        sides = 2, # 2022-07-21 00:20 -07:00 "side" in the original
        circular = TRUE
    )
  # reunite with pix
  splines = cbind(pix,splines)
  # returns a matrix, without the original column names
  # 2022-07-21 14:06 -07:00 but do we care? since it is 
  # called within the function to derive the smooth 
  # spline, and that returns scalars; not sure about
  # the derivatives, though
  return(splines)
} 

# filter is circular so "take" first and last values from 
# scaled data
# 2022-07-21 14:45 -07:00 by "take" is assumed to mean
# trim/remove
trim_scaled <- function(x) {
  # remove pix column, first and last values of the
  # scaled series
  trimmed = x[,-c(1,2,781)]
  return(trimmed)
}

# smooth spline on filtered data
# 2022-07-21 14:31 -07:00
# JDAY.x is loaded in global environment by the 
# prepare.R script via cons.R, which brings in obj/JDAY.x.rds
# 2022-07-21 14:58 -07:00 JDAY.x and x have differing
# lengths, due to removing first and last values in the
# smoothed series: 780 vs 778, so trim first  and last value of JDAY.x
# JDAY.x = JDAY.x[-1]
# x is a data.frame, so needs to be unlisted
# adjust for pix column by subsetting out column 1
# remove first and last values per lines 54-55 of
# dopixel_3g_v2_01.R
# pre-allocate a matrix sized to hold a result the same size as
# as the scaled object being provided as x
make_smooth <- function(x){
  # adjust col width to reflect removal of first and last values
  ss <- matrix(0,chunk_size,col_width - 2)
  for(i in 1:chunk_size)  {
      # TODO fix hardwired 781
      interior = x[i,-c(1:2,781)]
      smoothed[i,] = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                                   interior)$y
  }
  # Error in smoothed[i, ] <- smooth.spline(JDAY.x[-c(2, col_width - 1)],  :   # number of items to replace is not a multiple of replacement length
  return(ss)
}

# does work
# smooth.spline(JDAY.x[-c(2,col_width-1)],interior)$y


# WIP
#   lq <- quantile(ss$y, 0.1)
#   uq <- quantile(ss$y, 0.9)
#   amplitude <- uq - lq
#   # the mean evi of the time series
#   mean.evi <- mean(ss$y) 
#   # the standard deviation of the time series
#   sd.evi <- sd(ss$y) 
#   # sum evi
#   sum.evi <- sum(ss$y) 
#   
# # smoother spline and derivatives
# # 2022-07-21 22:26 -07:00
# # use predict.smooth.spline for the derivatives?
# 
#   sss <- smooth.spline(JDAY.x, approx.e, spar = 0.3)
#   # 2022-07-21 14:16 -07:00
#   # these support on.day from which the elon objects are derived
#   d1 <- predict(sss, sss$x, deriv = 1)$y
#   d1 <- (d1 - mean(d1)) / sd(d1)
#   d2 <- predict(sss, sss$x, deriv = 2)$y
#   d2 <- (d2 - mean(d2)) / sd(d2)
# 
#   devi.ss<-mean(abs(approx.e-ss$y)) # deviation between less smooth spline and data
#   devi.sss<-mean(abs(approx.e-sss$y))  # deviation between smoother spline and data
#   # need to cbind these with pix
# 
# # scale source data
# scale_ndvi <- function(x) {
#   # save out record id, pix
#   pix = x[,1]
#   x   = x[,-c(1)]
#   e = floor(x/10)/1000
#   r = x - floor(x/10)*10 +1
#   ifelse(e <= 0,NA,0)
#   e[which(r == 4 | r == 6)] = 0
#   e[which(r == 7)] = NA
#   # restore pix
#   e = cbind(pix,e)
#   return(e)
# }

#how many NA data points on scaled data needs used row-by-row
count_na <- function(x) length(which(is.na(x)))

# how many zero points on scaled data needs used row-by-row
count_zero <- function(x) length(which(x == 0))
