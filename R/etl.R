# etl.R
# extract transform load
# author: Richard Careaga
# Date: 2022-07-20

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# DATA

# outside data pre-processing

# install Rust https://www.rust-lang.org/ then see
# https://github.com/BurntSushi/xsv#installation for information
# from command line
# cargo install xsv # rust installation required
# xsv index d.csv
# xsv split splits d.csv

# fetch one first batch of 500 only 0.014 sec; less than a minute for total 
# input time, therefore

tic()
o <- fread(here("obj/splits/1000.csv"))
toc()

# apply the scaling function, scale_ndvi

# 500x780 data.table of scaled data returned in 0.5 sec
# about 30 minutes to handle the whole thing
# might be parallel, but need record id, because currently
# the row id is the implicit id
# scale_ndvi(d)  |> str()
# tic() result r is a data.table/data.frame
r <- scale_ndvi(o)
# toc()

# we can write back quite quickly with
fwrite(r, file = "scaled.csv", append = TRUE)

# if we do this 4,000 times, probably not too long

# assume we are now working with reading back in
# as the scaled dataset, now called s, in chunks of 500

s <- fread(here::here("scaled.csv"))

# plan: 
# 1. design output object
# 2. write functions to write to output object
# 3. test on a few chunks

# peak.day<-day.numbers[by(sss$y,YEAR,which.max)]
# trough.day<-day.numbers[by(sss$y,YEAR,which.min)]
# peak.day<-round( circ.mean(peak.day[-1]))
# trough.day<-round( circ.mean(trough.day[-1]))
# approx.e is s
# approx.e <- approx(JDAY.x,data.e,xout=JDAY.x,rule=2)$y
# 
splines <- make_splines(s)

make_smooth <- function(x){
  # adjust col width to reflect removal of first and last values
  #ss = matrix(0,chunk_size,col_width - 2) 
  for(i in 1:chunk_size)  {
    # TODO fix hardwired 781
    ss = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                         x[i,-c(1:2,781)])$y
    # TODO: need to capture pix, not i
    pix = i
    lq = quantile(ss, 0.1)
    uq = quantile(ss, 0.9)
    amplitude = uq - lq
    # TODO: write to receiver object
    # the mean evi of the time series
    mean.evi = mean(ss)
    # the standard deviation of the time series
    sd.evi = sd(ss)
    # sum evi
    sum.evi = sum(ss)
    
    results = c(pix,lq,uq,amplitude,mean.evi,sd.evi,sum.evi)
    names(results) <- c("pix","lq","uq","amplitude","mean.evi","sd.evi","sum.evi")
    return(results)
  }
}



