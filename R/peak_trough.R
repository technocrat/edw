# peak_trough.R
# develop function to calculate peak/trough
# author: Richard Careaga
# Date: 2022-08-24

# SETUP

source(here::here("R/prepare.R"))

# DATA

# first 2000 rows of unscaled data
# creates data.table dim 2000 781
d <- fread(here("obj/chunks/chunk_1.csv"))

# PREPROCESSING

# remove 'pix` identifier column
d <- d[,-1]
# use base matrix, as data.table::as.matrix transposes
d <- base::as.matrix(d,)

# MAIN

# Scale raw ndvi data by factor of 1e-4

# scale_ndvi converts to 780 2000, so transpose
d <- t(apply(d,1,scale_ndvi))

# change 5e5 NA values to doubles
# length(which(is.na(d)))
# returns equivalent of approx.e, but mean ex-NA is NaN
# transpose to 2000 780
d <- t(apply(d,1,interpolate))

# save out a copy for use in calculating sss
# d will be further pre-processed and used to
# calculate ss
approx.e <- d

# used after application of filter() in the
# make_splines() function, which is circular
begin = d[,1]
finish = d[,col_width]

# obtain filtered running mean
d <- make_splines(d)

# filter, used in make_splines, is circular, so replace
# first and last values by begin and finish
d[,1] <- begin
d[,col_width] <- finish

# smooth spline on filtered data; used to derive
# lq up amp mean.evi st.evi sum.evi

ss = smooth.spline(JDAY.x, t(d[1,]), all.knots = TRUE)$y

# deviation between less smooth spline and data
devi.ss <- apply(approx.e,1,find_mean_abs)

# smoother spline
sss_y <- t(apply(d,1,get_sss_y))

# deviation between smoother spline and data
devi.sss <- apply(approx.e,1,find_mean_abs)

# calculate peak.day, trough.day

# day.numbers <- sort(unique(JDAY))
# JDAY from cons.R can't be used because
# those are the unique ydays for all years of
# all pixels, not for any pixel or a year in it
# the_dates object, loaded from cons.R has the
# 780 dates in the series that apply to each pixel

# TODO:  how to apply this logic to object
# the_dates is a length 780 vector of Date objects,
# which can't co-exist in a matrix with the double

# objects in sss_y
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

trough.x <- ifelse(trough.day > 188, trough.day, trough.day + 365)
trough.x <- as.matrix(seq(trough.x, by = 365, length.out = 33))
window <- max(16, round((180*f[1,amp])/2, 0))
trough.win.x <- apply(trough.x, 1, FUN=function(x) (x-window):(x+window))
if(typeof(trough.win.x)=="integer") trough.win.x = rbind(trough.win.x)

td <- vector(length = 33)
td.x <- vector(length = 33)
td.evi <- vector(length = 33)

for (i in 1:33) {
  td.win = which(JDAY.x %in% trough.win.x[, i])
  if (length(td.win) > 0) {
    trough = td.win[which.min(sss$y[td.win])]
    td[i] = JDAY[trough]
    td.x[i] = JDAY.x[trough]
    td.evi[i] = ss[trough]
  } else {
    td[i] = NA
    td.x[i] = NA
    td.evi[i] = NA
  }
}

# needs to be in loop 
pheno.yr.win <- which(JDAY.x %in% (td.x[i]:(td.x[i+1])-1))
