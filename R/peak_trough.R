# peak_trough.R
# develop function to calculate peak/trough
# author: Richard Careaga
# Date: 2022-08-24

# TODO:  Identify code that takes the dim 2000 780 input chunks 
# as its argument and the code that requires the length 780 vector
# ie., matrix vs vector

# SETUP

source(here::here("R/prepare.R"))

# DATA

# first 2000 rows (pixels) of unscaled data
# creates data.table dim 2000 781
# I use d to signify data imported from source
# MATRIX
d <- fread(here("obj/chunks/chunk_1.csv"))

# PREPROCESSING

# remove 'pix` identifier column
# MATRIX
d <- d[,-1]

# MAIN

## Scale raw ndvi data by factor of 1e-4

# scale_ndvi converts to 780 2000, so transpose


d <- d[, lapply(.SD, scale_ndvi)]

# interpolate returns the equivalent of approx.e

## interpolate missing values

# MATRIX
d <- apply(d,1,interpolate)

# interpolate() transposes, so restore
d <- t(d)

# save out a copy for use in calculating sss
# d will be further pre-processed and used to
# calculate ss

# MATRIX 
approx.e <- d

# used after application of filter() in the
# make_splines() function, which is circular
# retain first and last values of each pixel
# MATRIX
begin = d[,1]
finish = d[,col_width]

##  obtain filtered running mean
# MATRIX
d <- make_splines(d)

# filter, used in make_splines, is circular, so replace
# first and last values by begin and finish
# MATRIX
d[,1] <- begin
d[,col_width] <- finish

## smooth spline on filtered data
# MATRIX but transposes
ss <- apply(d,1,smooth_spline)
ss <- t(ss)

# convert d back to data.table to add variables
# MATRIX
d <- data.table(d)

# TODO:  Verify use of quantiles rather than max/min
# MATRIX
d[,lq := quantile(ss, 0.1)]
d[,uq := quantile(ss, 0.9)]
d[,amplitude := uq - lq]
d[,mean.evi := mean(ss)]
d[,sd.evi := sd(ss)]
d[,sum.evi := sum(ss)]

# deviation between less smooth spline and data
# MATRIX returns length 2000 VECTOR
devi.ss <- apply(approx.e,1,find_mean_abs_ss)

# deviation between smoother spline and data
# MATRIX returns length 2000 VECTOR
devi.sss <- apply(approx.e,1,find_mean_abs_sss_y)

# MATRIX transposes to 780 2000, on purpose
sss_y <- apply(approx.e,1,get_sss_y)

# # save out copy for use in peak/trough
sss_Y <- data.table(sss_y)

# convert to data.table because dates and doubles in the 
# same matrix will convert to typeof character
# TODO:  Verify use of quantiles rather than max/min

# MATRIX
# TODO:  check doubled date columns
sss_y <- data.table(sss_y)
sss_y <- cbind(Date = ymd(the_dates),sss_y)
sss_y[, jday:=yday(Date)]
sss_y[, Year:=year(Date)]
sss_y[,year_min := min(jday),Year]
# TODO:  remind myself why we are interested in which
# day is the minimum jday for the year
the_min_rows <- which(sss_y$jday == sss_y$year_min)
sss_y[the_min_rows,mean(jday)]
sss_y[,year_max := max(jday),Year]
the_max_rows <- which(sss_y$jday == sss_y$year_max)
sss_y[the_max_rows,mean(jday)]

## calculate peak.day, trough.day

# MATRIX 

sss_Y[,Year := YEAR]
sss_Y[,Date := the_dates]
sss_Y[,jday := yday(the_dates)]
sss_Y <- sss_Y[,c(2001:2003,1:2000)]

# find the maximum by year and jday, 
peak.trough <- sss_Y[,lapply(.SD,max),by = c("Year","jday"), .SDcols = !c("Date")]

# select the jday of the maximum/minimum for each year
peak.day   <- peak.trough[,lapply(.SD,max),by = Year]
trough.day <- peak.trough[,lapply(.SD,min),by = Year] 

# TODO:  Understand what this should be doing;

peak.day$year_max <- round(apply(peak.day[,-1],1,circ.mean))
trough.day$year_min <- round(apply(trough.day[,-1],1,circ.mean))

# smoother spline and derivatives

# approx.e is t() transposed because it is dim 780 2000
# which means that each row has length 2000, while JDAY.x 
# used in make_sss has length 780
# ds contains lists for d1 and d2, each of which is a list of lists
ds <- apply(t(approx.e),1,make_sss)

trough.x <-  adjust_trough.day(trough.day$year_min)

# subset of trough.x is necessary because trough.x is a vector,of 
# length 33 and the first argument to seq must be length 1
trough.x <- as.matrix(seq(trough.x[1], by = 365, length.out = 33))

# window is determined for each pixel individually per Higgens
# p. 3584, column 1

# changed from window, which is a closure
win <- max(16,round(d[,amplitude] * 180 / 2,0))

# TODO:  Fix--this won't work because FUN creates a dimensional-less
# vector, and apply expects a dimension
trough.win.x <- apply(trough.x, 1, FUN=function(x) (x-win):(x+win))

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
