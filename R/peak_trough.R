# peak_trough.R
# develop function to calculate peak/trough
# author: Richard Careaga
# Date: 2022-08-24

# SETUP

source(here::here("R/prepare.R"))

# DATA

# first 2000 rows of unscaled data
# creates data.table dim 2000 781
# I use d to signify data imported from source
d <- fread(here("obj/chunks/chunk_1.csv"))

# PREPROCESSING

# remove 'pix` identifier column
d <- d[,-1]


# use base matrix, as data.table::as.matrix transposes
#d <- base::as.matrix(d,)

# MAIN

## Scale raw ndvi data by factor of 1e-4

# scale_ndvi converts to 780 2000, so transpose

d <- d[, lapply(.SD, scale_ndvi)]

# interpolate returns the equivalent of approx.e

## interpolate missing values

d <- apply(d,1,interpolate)
# save out a copy for use in calculating sss
# d will be further pre-processed and used to
# calculate ss

approx.e <- d

# used after application of filter() in the
# make_splines() function, which is circular
begin = d[,1]
finish = d[,col_width]

##  obtain filtered running mean
d <- make_splines(d)

# filter, used in make_splines, is circular, so replace
# first and last values by begin and finish
d[,1] <- begin
d[,col_width] <- finish

## smooth spline on filtered data
ss <- apply(t(d),1,smooth_spline)

# convert d back to data.table to add variables
d <- data.table(d)
# TODO:  Verify use of quantiles rather than max/min

d[,lq := quantile(ss, 0.1)]
d[,uq := quantile(ss, 0.9)]
d[,amplitude := uq - lq]
d[,mean.evi := mean(ss)]
d[,sd.evi := sd(ss)]
d[,sum.evi := sum(ss)]

# deviation between less smooth spline and data
devi.ss <- apply(approx.e,1,find_mean_abs_ss)

# deviation between smoother spline and data
devi.sss <- apply(t(approx.e),1,find_mean_abs_sss_y)

sss_y <- apply(t(approx.e),1,get_sss_y)

# convert to data.table because dates and doubles in the 
# same matrix will convert to typeof character
# TODO:  Verify use of quantiles rather than max/min

sss_y <- data.table(sss_y)
sss_y <- cbind(Date = ymd(the_dates),sss_y)
sss_y[, jday:=yday(Date)]
sss_y[, Year:=year(Date)]
sss_y[,year_min := min(jday),Year]
the_min_rows <- which(sss_y$jday == sss_y$year_min)
sss_y[the_min_rows,mean(jday)]
sss_y[,year_max := max(jday),Year]
the_max_rows <- which(sss_y$jday == sss_y$year_max)
sss_y[the_max_rows,mean(jday)]

## calculate peak.day, trough.day

# TODO:  Why is there a single time series for all pixels?

peak.day <- unique(sss_y[,year_max, by = "YEAR"])
trough.day <- unique(sss_y[,year_min, by = "YEAR"])

peak.day <- round(apply(peak.day[,-1],1,circ.mean))
trough.day <- round(apply(trough.day[,-1],1,circ.mean))

# smoother spline and derivatives

# ds contains lists for d1 and d2, each of which is a list of lists
ds <- apply(approx.e,1,make_sss)

trough.day <-  adjust_trough.day(trough.day$year_min)

trough.x <- as.matrix(seq(trough.x[1], by = 365, length.out = 33))

window <- max(16, round((180*amplitude)/2, 0))
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
