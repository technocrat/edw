# scratch.R
# scratch pad
# author: Richard Careaga
# Date: 2022-08-11

source(here::here("R/prepare.R"))

# first 2000 rows of unscaled data
d <- fread(here("obj/chunks/chunk_1.csv"))
# remove 'pix` identifier column
d <- d[, -1]
d <- as.matrix(d)
d <- scale_ndvi(d)

d <- t(apply(d, 1, interpolate))
approx.e <- d
begin <- d[, 1]
finish <- d[, 780]
d <- make_splines(d) # d is equivalent of data.e
data.e <- d

get_sss_x <- function(x) smooth.spline(JDAY.x, x, spar = 0.3)$x
get_sss_y <- function(x) smooth.spline(JDAY.x, x, spar = 0.3)$y

# dim 780 2000--ie dates are rows and columns are pixels
a <- apply(data.e, 1, get_sss_y)

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

# YEAR from cons.R
# by is a convenience wrapper for tapply
# sss.y is a 780 2000 matrix typeof double
#
# YEAR should be a factor?
# peak.day <- day.numbers[by(sss$y, YEAR, which.max)]
# trough.day <- day.numbers[by(sss$y, YEAR, which.min)]
# peak.day <- round(circ.mean(peak.day[-1]))
# trough.day <- round(circ.mean(trough.day[-1]))

# all years, pixel 250 of chunk
b <- a[, 250]
# can't be matrix--b/c converts to character
f <- cbind(the_dates, b)
f <- data.table(the_dates, b)
f[, jday := yday(the_dates)]
f[, Year := year(the_dates)]
day.numbers <- unique(f$jday)
f[, year_min := min(b), Year]
the_min_rows <- which(f$b == f$year_min)
# this gives DT w/33 rows, with the jday for each year of
# one pixel
f[the_min_rows, ]

mean_trough <- f[, mean(year_min)]
f[, year_max := max(b), Year]
# the_max_rows <- which(f$b == f$year_max)
# # this gives DT w/33 rows, with the jday for each year of
# # one pixel
# # arithmetic
# f[the_max_rows,mean(jday)]
# # 233.1818
# # library(circular) no difference
# f[the_max_rows,mean(jday)]
f[f[, b == year_min]]

# f has already been splined, but not smooth splined

# filter is circular so take first value from data
begin <- f[1, 1]
finish <- f[780, 1]
# f used for find_mean_abs and make_sss
# fs is used for the smooth.spline derived objects,
# lq uq amplitude, mean.evi,sd.evi,sum.evi
# fs is ss in dopixel
fs <- smooth.spline(JDAY.x, f[, b], all.knots = TRUE)$y
lq <- quantile(fs, 0.1)
uq <- quantile(fs, 0.9)
amplitude <- uq - lq
mean.evi <- mean(fs)
sd.evi <- sd(fs)
sum.evi <- sum(fs)

# approx e is return value of interpolate on chunk 2000 780
approx.e <- approx.e[1, ]
# smoother spline and deriviatives
sss <- smooth.spline(x = JDAY.x, y = t(approx.e), spar = 0.3)
d1 <- predict(sss, sss$x, deriv = 1)$y
d1 <- (d1 - mean(d1)) / sd(d1)
d2 <- predict(sss, sss$x, deriv = 2)$y
d2 <- (d2 - mean(d2)) / sd(d2)

peak.day <- day.numbers[by(sss$y, YEAR, which.max)]
trough.day <- day.numbers[by(sss$y, YEAR, which.min)]
peak.day <- round(circ.mean(peak.day[-1]))
trough.day <- round(circ.mean(trough.day[-1]))

trough.x <- ifelse(trough.day > 188, trough.day, trough.day + 365)
trough.x <- as.matrix(seq(trough.x, by = 365, length.out = 33))

window <- max(16, round((180 * amplitude) / 2, 0))
trough.win.x <- apply(trough.x, 1, FUN = function(x) (x - window):(x + window))
# error: "condition has length > 1"
# if(class(trough.win.x)=="integer") trough.win.x<-rbind(trough.win.x)
if (typeof(trough.win.x) == "integer") trough.win.x <- rbind(trough.win.x)

td <- vector(length = 33)
td.x <- vector(length = 33)
td.evi <- vector(length = 33)
for (i in 1:33) {
  td.win = which(JDAY.x %in% trough.win.x[, i])
  if (length(td.win) > 0) {
    trough = td.win[which.min(sss$y[td.win])]
    td[i] = JDAY[trough]
    td.x[i] = JDAY.x[trough]
    td.evi[i] = fs[trough]
  } else {
    td[i] = NA
    td.x[i] = NA
    td.evi[i] = NA
  }
}
