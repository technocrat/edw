# peak_trough.R
# develop function to calculate peak/trough
# author: Richard Careaga
# Date: 2022-08-24

# SETUP

source(here::here("R/prepare.R"))

# DATA

# first 2000 rows (pixels) of unscaled data
# creates data.table dim 2000 781
# I use d to signify data imported from source

d <- fread(here("obj/chunks/chunk_1.csv"))

# PREPROCESSING

# remove 'pix` identifier column

d <- d[,-1]

# MAIN

## Scale raw ndvi data by factor of 1e-4

d <- d[, lapply(.SD, scale_ndvi)]

## interpolate missing values
# interpolate returns the equivalent of approx.e

d <- apply(d,1,interpolate)

# save out a copy for use in calculating sss
# d will be further pre-processed and used to
# calculate ss

approx.e <- d

# used after application of filter() in the
# make_splines() function, which is circular
# retain first and last values of each pixel

begin = d[,1]
finish = d[,col_width]

##  obtain filtered running mean

d <- make_splines(d)

# filter, used in make_splines, is circular, so replace
# first and last values by begin and finish

d[,1] <- begin
d[,col_width] <- finish

## smooth spline on filtered data
# TODO:  Should these be so small?

# transpose to have rows as pixels
d <- t(d)
ss <- apply(d,1,smooth_spline)

# convert d back to data.table to add variables
d <- data.table(d)

# TODO:  Verify use of quantiles rather than max/min

d[,lq := apply(ss,2,q1)]
d[,uq := apply(ss,2,q9)]
d[,amplitude := uq - lq]
d[,mean.evi := apply(ss,2,mean)]
d[,sd.evi := apply(ss,2,sd)]
d[,sum.evi := apply(ss,2,sum)]

# deviation between less smooth spline and data

devi.ss <- apply(approx.e,1,find_mean_abs_ss)

# deviation between smoother spline and data

devi.sss <- apply(approx.e,2,find_mean_abs_sss_y)


sss_y <- apply(approx.e,2,get_sss_y)

# save out copy for use in peak/trough
sss_Y <- data.table(sss_y)

# convert to data.table because dates and doubles in the 
# same matrix will convert to typeof character

# TODO:  Verify use of quantiles rather than max/min

sss_y <- data.table(sss_y)
sss_y <- cbind(Date = ymd(the_dates),sss_y)
sss_y[, jday:=yday(Date)]
sss_y[, Year:=year(Date)]
sss_y[,year_min := min(jday),Year]
the_min_rows <- which(sss_y$jday == sss_y$year_min)
sss_y[,year_max := max(jday),Year]
the_max_rows <- which(sss_y$jday == sss_y$year_max)

## calculate peak.day, trough.day

sss_Y <- sss_y
sss_Y[,Year := YEAR]
sss_Y[,Date := the_dates]
sss_Y[,jday := yday(the_dates)]
sss_Y <- sss_Y[,c(2001:2003,1:2000)]

# find the maximum by year and jday, 

peak.trough <- sss_Y[,lapply(.SD,max),by = c("Year","jday"), .SDcols = !c("Date")]

# select the jday of the maximum/minimum for each year

peak.day   <- peak.trough[,lapply(.SD,max),by = Year]
trough.day <- peak.trough[,lapply(.SD,min),by = Year] 

peak.day$year_max <- round(apply(peak.day[,-1],1,circ.mean))
trough.day$year_min <- round(apply(trough.day[,-1],1,circ.mean))
# smoother spline and derivatives
# ds contains lists for d1 and d2, each of which is a list of lists

ds <- apply(approx.e,2,make_sss)

trough.x <-  adjust_trough.day(trough.day$year_min)

# subset of trough.x is necessary because trough.x is a vector,of 
# length 33 and the first argument to seq must be length 1

trough.x <- as.matrix(seq(trough.x[1], by = 365, length.out = 33))

# window is determined for each pixel individually per Higgens
# p. 3584, column 1; so it needs to have a value for each pixel,
# not all pixels
# changed from `window`, which is a closure

adjusted_amplitudes <- round(d[,amplitude] * 180 / 2,0)
win <- ifelse(adjusted_amplitudes < 16,16,adjusted_amplitudes)
# 33 x 33 matrix 349:1208, column-wise from trough.x, which
# is 1x1 365:12045, same as `win`
# 
trough.win.x <- apply(trough.x, 1, FUN=function(x) (x-win):(x+win))

# TODO:   figure out what this is supposed to do, even though it evaluates to TRUE, only binds to itself
if(typeof(trough.win.x)=="integer") trough.win.x = rbind(trough.win.x)

td <- vector(mode = "numeric", length = 33)
td.x <- vector(mode = "numeric",length = 33)
td.evi <- vector(mode = "numeric",length = 33)

# this returns a list of 2000 lists of smooth.spline
# objects, which are lists of 21 objects
sss_list  <- apply(approx.e,2,get_sss)

for (i in 1:33) {
  td.win = which(JDAY.x %in% trough.win.x[, i])
  # 130 JDAY.x are there
  if (length(td.win) > 0) {
    # TODO:  why does this evaluate to 0?
    # b/c causes td[i] <- JDAY[trough] to fail
    # so need to fix
    # sss$y is NOT the same as sss_y
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
