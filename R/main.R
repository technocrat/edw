# main.R
# accumulate operations of csv files 
# author: Richard Careaga
# Date: 2022-08-03

# load libraries, functions and constants

source(here::here("R/prepare.R"))

# RECEIVER FOR FIRST BATCH

M = matrix(nrow = 386, ncol = 12)

# DATA

## FIRST BATCH

leftover <- as.matrix(fread(here("obj/chunks/chunk_1047.csv"))) # only 386 rows

M[,1] <- leftover[,1] # save pix
leftover <- leftover[,-1] # trim pix
chunk <- scale_ndvi(leftover)
# populate M
source(here("R/etl.R"))
# save to data.table
saveRDS(rbind(receiver,M, use.names = FALSE), file = here("obj/processed_leftover.Rds"))
# take leftover objects out of memory
rm(leftover,chunk,chunk_ss,M)

# REMAINING BATCHES

# files to iterate over

# each is 2000 rows, except for csv_1047.csv, which is 386
chunks <- dir(here("obj/chunks"), pattern = "\\.csv$", full.names = TRUE) 

# remove stub
forepart = paste0(here(),"/obj/chunks/chunk_")
aftpart <- paste0("1047.csv")
target <- paste0(forepart,aftpart)
chunks <- chunks[-which(chunks == target)]

# TESTING
# set.seed(42)
# # Took 25 sec, same 0.00125
# pieces <- sample(chunks,10,replace = FALSE)
# # Took 253 sec, same 00125
# pieces <- sample(chunks,100,replace = FALSE)


for(chunk in chunks) {
  M = matrix(nrow = 2000, ncol = 12)
  chunk = as.matrix(fread(here(chunks[1])))
  M[,1] = chunk[,1] # save pix
  chunk = chunk[,-1] # trim pix
  chunk = apply(d,1,scale_ndvi)

  # populate M
  source(here("R/etl.R"))
  receiver = rbind(receiver,M, use.names = FALSE)
}

leftover <- readRDS(here("obj/processed_leftover.Rds"))
receiver = rbind(receiver,leftover)
saveRDS(receiver, file = here("obj/processed_receiver.Rds"))

