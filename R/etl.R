# etl.R
# extract transform load
# author: Richard Careaga
# Date: 2022-07-20

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# DATA

chunks <- dir(here("obj/chunks"), pattern = "\\.csv$", full.names = TRUE)
leftover <- fread(here("chunks/chunk_1047.csv")) # only 387 lines

# ACCUMULATION OBJECTS

# brings empty data.table to take output of for loops below

receiver <- source(here("R/receiver.R"))
passer <- receiver
# exclude last chunk, which has an off number of rows

for(chunk in seq_along(chunks-1)) {
  for(i in seq_along(chunk)) {
    # accumulate results in a vector
    v = vector(length = 12)
    v[1] = chunk[,1] # pix value for this chunk
    the_chunk = scale_ndvi(chunk)
    # discard pix
    the_chunk = the_chunk[,-1]
    v[11] = length(which(is.na(the_chunk)))
    v[12] = length(
      which(the_chunk - 
              floor(the_chunk/10)*10 + 1 == 4 | 
              the_chunk - 
              floor(the_chunk/10)*10 + 1 == 6))
    the_chunk = interpolate(the_chunk)
    the_chunk = make_splines(the_chunk)
    v[2:8] = make_smooth(the_chunk)
    v[9] = make_sss(the_chunk)
    v[10] = L.JDAY
    # accumulate into data.table
    passer = passer
  }
}

# last chunk
for(chunk in seq_along(leftover)) {
  for(i in seq_along(chunk)) {
    
    the_chunk = scale_ndvi(chunk)
    pix_list = vector(mode = "list",length = chunk_size)
    pix_list[i] = the_chunk[i,1]
    the_chunk = the_chunk[,-1]
    the_na = length(which(is.na(the_chunk)))
    the_ice = length(
      which(the_chunk - 
              floor(the_chunk/10)*10 + 1 == 4 | 
              the_chunk - 
              floor(the_chunk/10)*10 + 1 == 6))
    the_chunk = interpolate(the_chunk)
    the_chunk = make_splines(the_chunk)
    the_smooth = make_smooth(the_chunk)
    make_sss(the_chunk)
  }
}
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

for(i in seq_along(v)) passer[1,i] = v[i]