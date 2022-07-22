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
# approx.e is 
# approx.e<-approx(JDAY.x,data.e,xout=JDAY.x,rule=2)$y
# 
splines <- make_splines(s)

sss <- smooth.spline(JDAY.x, x, spar=0.3)
