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

d <- t(apply(d,1,interpolate))
d <- make_splines(d)

get_sss_y <- function(x) smooth.spline(JDAY.x,x,spar = 0.3)$y
  
a <- apply(d,1,get_sss_y)


# calculate peak.day, trough.day

# JDAY from cons.R
day.numbers <- sort(unique(JDAY))

# YEAR from cons.R
# by is a convenience wrapper for tapply
# sss.y is a 780 2000 matrix typeof double
# 
# YEAR should be a factor?
# peak.day <- day.numbers[by(sss$y, YEAR, which.max)]
# trough.day <- day.numbers[by(sss$y, YEAR, which.min)]
# peak.day <- round(circ.mean(peak.day[-1]))
# trough.day <- round(circ.mean(trough.day[-1]))

# get first year
py1 <- a[1:13]
min(py1)
max(py1)
# second year is better example, because it is full year
py2 <- a[1,14:37]
py2.peak.day <- max(py2)
py2.trough.day <- min(py2)
for(y in YEAR) 
  
b <- a[,1]
# can't be matrix--b/c converts to character
f <- cbind(the_dates,b)
f <- data.table(the_dates,b)
f[, jday:=yday(the_dates)]
f[, Year:=year(the_dates)]
g <- f
# b is sss$y
g[,year_min := min(b),Year]
the_min_rows <- which(g$b == g$year_min)
# this gives DT w/33 rows, with the jday for each year of
# one pixel
g[the_min_rows,]

