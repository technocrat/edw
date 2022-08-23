# scratch.R
# scratch pad
# author: Richard Careaga
# Date: 2022-08-11

source(here::here("R/prepare.R"))

# first 2000 rows of unscaled data
d <- fread(here("obj/chunks/chunk_1.csv"))
# remove 'pix` identifier column
d <- d[,-1]
d <- as.matrix(d)
d <- scale_ndvi(d)

# returns equivalent of approx.e
d <- t(apply(d,1,interpolate))
approx.e <- t(d)
begin = d[,1]
finish = d[,780]
d <- make_splines(d)
# filter, used in make_splines, is circular, so replace
# first and last values by begin and finish
d[,1] <- begin
d[,col_width] <- finish
# d is now used for find_mean_abs and make_sss

# dim 780 2000--ie dates are rows and columns are pixels
sss_y <- apply(d,1,get_sss_y)

# calculate peak.day, trough.day
# First, the average Julian day of
# minimum NDVI (trough day) was calculated for each pixel by
# averaging the trough days of all years in the 31-year time ser-
#  ies
# NB: we have 32 years plus 1981, a half year, not 31

# JDAY from cons.R
# however, those are the unique ydays for all years of
# all pixels
# day.numbers <- sort(unique(JDAY))


# all years, pixel 250 of chunk
b <- sss_y[,250]
# can't be matrix--b/c converts to character
f <- cbind(the_dates,b)
f <- data.table(the_dates,b)
f[, jday:=yday(the_dates)]
f[, Year:=year(the_dates)]
day.numbers <- unique(f$jday)
f[,year_min := min(b),Year]
the_min_rows <- which(f$b == f$year_min)
# this gives DT w/33 rows, with the jday for each year of
# one pixel
f[the_min_rows,] 
f[the_min_rows,mean(jday)]
f[,year_max := max(b),Year]
the_max_rows <- which(f$b == f$year_max)
# this gives DT w/33 rows, with the jday for each year of
# one pixel
# arithmetic
f[the_max_rows,mean(jday)] 
# 233.1818
# library(circular) no difference
f[the_max_rows,mean(jday)]

#f has already been smooth splined, so we can find 
#uq lq amplitude directly

f[,lq := quantile(b, 0.1),Year]
f[,uq := quantile(b, 0.9),Year]
f[,amp := uq - lq,Year]

# YEAR from cons.R
# by is a convenience wrapper for tapply
# sss.y is a 780 2000 matrix typeof double

peak.day <- day.numbers[by(f$b, YEAR, which.max)]
trough.day <- day.numbers[by(f$b, YEAR, which.min)]
peak.day <- round(circ.mean(peak.day[-1]))
trough.day <- round(circ.mean(trough.day[-1]))

sss <- smooth.spline(JDAY.x, t(approx.e)[1,], spar = 0.3)
d1 <- predict(sss, sss$x, deriv = 1)$y
d1 <- (d1 - mean(d1)) / sd(d1)
d2 <- predict(sss, sss$x, deriv = 2)$y
d2 <- (d2 - mean(d2)) / sd(d2)