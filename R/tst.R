# tst.R
# test output of my code against output of provided code
# author: Richard Careaga
# Date: 2022-08-07

source(here::here("R/prepare.R"))

# first 2000 rows of unscaled data
d <- fread(here("obj/chunks/chunk_1.csv"))

# take first 500, removing pix
d <- as.matrix(d[1:500,-1])

# results provided
res0 <- readRDS(here("data/ndvi_df.Rds"))

# results produced
res1 <- readRDS(here("obj/processed_receiver.Rds"))
res1 <- res1[1:500,-1]
res1 <- as.matrix(res1)

# trim res0 to match column size of res1, which does not
# calculate as many variables
res0 <- res0[,c(1:6,9:13)]

# half the data are NA
# all non-NA L.JDAY are 780
summary(res0)
# no NA values, L.JDAY are all 7809
summary(res1)

# remove L.JDAY--no difference
res0 <- res0[,-9]
res1 <- res1[,-9]
  
# are there any NA in the data? No. Then why are there in res0?
   
DescTools::CountCompCases(d)

res0[,2] + res0[,1] == res0[,6]
res1[,2] + res1[,1] == res1[,6]

res0[1:50,c(2,1,6)]
res1[1:50,c(2,1,6)]
