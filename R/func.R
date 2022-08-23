# func.R
# functions used in scripts in this directory
# author: Richard Careaga
# Date: 2022-v07-20

circ.mean <- function(x) {
  x = rad(x * (360 / 365))
  sinr = sum(sin(x))
  cosr = sum(cos(x))
  circmean = atan2(sinr, cosr)
  circmean = deg(circmean)
  if (circmean < 0) circmean = 365 + circmean
  return(circmean)
}

rad <- function(degree) (degree * pi)/180

deg <- function(radian) (radian * 180)/pi

find_mean_abs <- function(x) mean(abs(x[-c(1,col_width - 1)] - ss))

get_quantiles <- function(x) quantile(x,probs = c(0.1,0.9),na.rm = TRUE)

# interpolate missing values of scaled data
# takes a single row of the scaled data


interpolate <- function(x) t(approx(JDAY.x,x, xout = JDAY.x, rule = 2)$y)

# takes a single row of the scaled data

make_splines <- function(x) {
  # get a filtered running mean: 2022-07-30 14:17 -07:00
  # transpose is to make nrow(x) > 1, since we process one
  # row at a time and the length of filter must be greater
  # than the number of rows
  # produces data.e in dopixel script
  splines = filter(t(x),
        filter = rep(1 / 3,3),
        method = "convolution",
        sides = 2, # 2022-07-21 00:20 -07:00 "side" in the original
        circular = TRUE
    )
  return(t(splines))
} 

# takes a single row of the scaled data

make_smooth <- function(x){
  #smoothed spline on filtered data
  ss = x
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

scale_ndvi <- function(x) {
  e = floor(x/10)/1000
  r = x - floor(x/10)*10 +1
  e = ifelse(e <= 0,NA,e)
  e[which(r == 4 | r == 6)] = 0
  e[which(r == 7)] = NA
  return(e)
}

get_sss_both <- function(x) smooth.spline(JDAY.x,x,spar = 0.3)
  