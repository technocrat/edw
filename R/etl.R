# etl.R
# extract transform load
# author: Richard Careaga
# Date: 2022-07-20

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# DATA

chunks <- dir(here("obj/chunks"), pattern = "\\.csv$", full.names = TRUE)
leftover <- fread(here("obj/chunks/chunk_1047.csv")) # only 387 lines

# ACCUMULATION OBJECTS

# brings empty data.table to take output of for loops below

receiver <- source(here("R/receiver.R"))

passer <- data.table(
  # pix is the record id added corresponding to original row order
  pix = 0,
  # variables derived from the smoothed splines
  lq = 0,0,
  uq = 0.0,
  mean.evi = 0.0,
  sd.evi = 0.0,
  sum.evi = 0.0,
  amplitude = 0.0,
  devi.ss = 0.0,
  devi.sss = 0.0,
  # a constant from cons.R
  L.JDAY = 0,
  # variables derived from from scaled data
  NA.length = 0,
  ICE.length = 0
)

 %>% 