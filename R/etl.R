# etl.R
# extract transform load
# author: Richard Careaga
# Date: 2022-08-03

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# RECEIVER FOR FIRST BATCH

M = matrix(nrow = 386, ncol = 12)

# DATA

## FIRST BATCH

leftover <- as.matrix(fread(here("obj/chunks/chunk_1047.csv"))) # only 387 rows

M[,1] <- leftover[,1] # save pix
leftover <- leftover[,-1] # trim pix
chunk <- scale_ndvi(leftover)
# populate M
source(here("R/main.R"))
# save to data.table
saveRDS(rbind(receiver,M, use.names = FALSE), file = here("obj/processed_leftover.Rds"))
# take leftover objects out of memory
rm(leftover,chunk,chunk_ss,M)

# REMAINING BATCHES

# files to iterate over

chunks <- dir(here("obj/chunks"), pattern = "\\.csv$", full.names = TRUE) # each is 2000 rows
library(tictoc)
tic()
#for(chunk in chunks[1]) {
  M = matrix(nrow = 2000, ncol = 12)
  chunk = as.matrix(fread(here(chunks[1])))
  M[,1] = chunk[,1] # save pix
  chunk = chunk[,-1] # trim pix
  chunk = scale_ndvi(chunk)
  # populate M
  source(here("R/main.R"))
  receiver = rbind(receiver,M, use.names = FALSE)
#}
toc()