# func.R
# functions used in scripts in this directory
# author: Richard Careaga
# Date: 2022-07-20

# interpolate missing values of scaled data
# takes a single row of the scaled data
 
interpolate <- function(x) approx(JDAY.x,x, xout = JDAY.x, rule = 2)$y

# takes a single row of the scaled data

make_splines <- function(x) {
  # get a filtered running mean:
  splines = filter(x,
        filter = rep(1 / 3,3),
        method = "convolution",
        sides = 2, # 2022-07-21 00:20 -07:00 "side" in the original
        circular = TRUE
    )
  return(splines)
} 

# takes a single row of the scaled data

make_smooth <- function(x){
  ss = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                     x[-c(1,col_width - 1)])$y
  lq = quantile(ss, 0.1)
  uq = quantile(ss, 0.9)
  amplitude = uq - lq
  # the mean evi of the time series
  mean.evi = mean(ss)
  # the standard deviation of the time series
  sd.evi = sd(ss)
  # sum evi
  sum.evi = sum(ss)
  # deviation between less smooth spline and data
  devi.ss <- mean(unlist(abs(x[-c(1,col_width - 1)]-ss))) 
  results = c(lq,uq,amplitude,mean.evi,sd.evi,sum.evi,devi.ss)
  names(results) <- c("lq","uq","amplitude","mean.evi","sd.evi","sum.evi","devi.ss")
  return(results)
}

# takes a single row of the scaled data

make_sss <- function(x){
  sss_both = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                           x[-c(2,col_width -1)],
                           spar = 0.3)
  sss.x = sss_both$x
  sss.y = sss_both$y
  d0 = predict(sss_both, sss.x, deriv = 1)$y
  d1 = (d0 - mean(d0)) / sd(d0)
  d0 = predict(sss_both, sss.x, deriv = 2)$y
  d2 = (d0 - mean(d0)) / sd(d0)
  # deviation between smoother spline and data
  devi.sss = mean(abs(x[-c(2,col_width -1)] - sss.y))
  return(devi.sss)
}

# save the current chunk's pix column

save_pix <- function(x) x[,1]

# scale source data (not row-by-row)

scale_ndvi <- function(x) {
  # save out record id, pix
  pix = x[,1]
  x   = x[,-c(1)]
  e = floor(x/10)/1000
  r = x - floor(x/10)*10 +1
  ifelse(e <= 0,NA,0)
  e[which(r == 4 | r == 6)] = 0
  e[which(r == 7)] = NA
  # restore pix
  e = cbind(pix,e)
  return(e)
}