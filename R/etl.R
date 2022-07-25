# etl.R
# extract transform load
# author: Richard Careaga
# Date: 2022-07-20

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# DATA

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

splines <- make_splines(s)





