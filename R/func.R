# func.R
# functions used in scripts in this directory
# author: Richard Careaga
# Date: 2022-08-29

adjust_trough.day <- function(x) ifelse(x > 188, x, x + 365)

circ.mean <- function(x) {
  x = rad(x * (360 / 365))
  sinr = sum(sin(x))
  cosr = sum(cos(x))
  circmean = atan2(sinr, cosr)
  circmean = deg(circmean)
  if (circmean < 0) circmean = 365 + circmean
  return(circmean)
}

find_r_value <- function(x) x - floor(x/10) * 10 + 1

find_mean_abs_ss <-    function(x) mean(abs(x - ss))
find_mean_abs_sss_y <- function(x) mean(abs(x - get_sss(x)$y))

get_sss <- function(x) smooth.spline(JDAY.x, x, spar=0.3)
get_sss_y <- function(x) smooth.spline(JDAY.x, x, spar=0.3)$y

get_quantiles <- function(x) quantile(x,probs = c(0.1,0.9),na.rm = TRUE)

# interpolate missing values of scaled data

interpolate <- function(x) approx(JDAY.x,x, xout = JDAY.x, rule = 2)$y

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

make_sss <- function(x){
  sss = get_sss(x)
  d1 = predict(sss, sss$x, deriv = 1)$y
  d1 = (d1 - mean(d1)) / sd(d1)
  d2 = predict(sss, sss$x, deriv = 2)$y
  d2 = (d2 - mean(d2))/sd(d2)
  return(list(d1,d2))
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