# prep.R
# convert data from ndvi.bin ndvi.desc to csv and save
# author: Richard Careaga
# Date: 2022-07-20

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# DATA

input <- bigmemory::attach.big.matrix(here("data/ndvi.desc"))

# PREPROCESSING

# extract to a data table object

d <- as.data.table(input[,], sorted = FALSE)

# MAIN

# add a record id based on pixel row number

pix <- as.data.table(c(pix = 1:size))

d <- cbind(pix,d)
colnames(d)[1] <- "pix"

# save d as a csv file

fwrite(d,here("obj/d.csv"))




