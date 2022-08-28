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

rad <- function(x) (x * pi)/180

deg <- function(x) (x * 180)/pi

find_r_value <- function(x) x - floor(x/10) * 10 + 1

find_mean_abs_ss <-    function(x) mean(abs(x - ss))
find_mean_abs_sss_y <- function(x) mean(abs(x - smooth.spline(JDAY.x,x,spar = 0.3)$y))

get_sss <- function(x) smooth.spline(JDAY.x, x, spar=0.3)

get_quantiles <- function(x) quantile(x,probs = c(0.1,0.9),na.rm = TRUE)

# interpolate missing values of scaled data

interpolate <- function(x) t(approx(JDAY.x,x, xout = JDAY.x, rule = 2)$y)

make_splines <- function(x) {
  # get a filtered running mean: 2022-07-30 14:17 -07:00
  # transpose is to make nrow(x) > 1, since we process one
  # row at a time and the length of filter must be greater
  # than the number of rows
  # produces data.e in dopixel script
  splines = filter(t(x),
        filter = rep(1 / 3,3),
        method = "convolution",
        # 2022-07-21 00:20 -07:00 "side" in the original
        sides = 2, 
        circular = TRUE
    )
  return(t(splines))
} 

make_smooth <- function(x){
  #ss is the smoothed spline on filtered data
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

make_sss <- function(x){
  sss.x = sss_both$x
  sss.y = sss_both$y
  d0 = predict(sss_both, sss.x, deriv = 1)$y
  d1 = (d0 - mean(d0)) / sd(d0)
  d0 = predict(sss_both, sss.x, deriv = 2)$y
}

scale_ndvi <- function(x) {
  e = floor(x/10)/1000
  r = x - floor(x/10)*10 +1
  e = ifelse(e <= 0,NA,e)
  e[which(r == 4 | r == 6)] = 0
  e[which(r == 7)] = NA
  return(e)
}

smooth_spline <- function(x) smooth.spline(JDAY.x,x,all.knots = TRUE)$y