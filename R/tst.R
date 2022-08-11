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
# calculate as peak.day trough. or any variables after
# ICE.length
res0 <- res0[,c(1:6,9:13)]

# half the provided results  are NA
# all non-NA L.JDAY are 780
summary(res0)
# no NA values, L.JDAY are all 780
summary(res1)

# remove L.JDAY--no difference
res0 <- res0[,-9]
res1 <- res1[,-9]
  
# are there any NA in the data? No. Then why are there in res0?
   

d_na <- DescTools::CountCompCases(d)
mean(d_na$tab$nas)

# amplitude = uq - lq 
f0 <- function(x) (x[2] - x[1]) == x[6]
apply(res0,1,f0)


# all wrong; TODO: revise calculation of amplitude for rowwise
apply(res1,1,f0)

# Does scaling introduce NAs?

# scale_ndvi operates on the raw data, here `d`

chunk <- scale_ndvi(d)

# however, dopixel operates rowwise

  # data.e<-floor(d[p,]/10)/1000
  # data.r<-data[p,]-floor(d[p,]/10)*10 +1
  # data.e[data.e<=0]<-NA
  # data.e[which(data.r==4 | data.r==6)]<-0
  # data.e[data.r==7]<-NA

# ∴ correct scale_ndvi() to operate rowwise
# create `chunk` scaled `d`
# TODO: why is chunk transposed?
chunk <- apply(d,1,scale_ndvi) # in place of scale_ndvi(d)
chunk <- t(chunk)

# OK, scaling does introduce NA
chunk_na <- DescTools::CountCompCases(chunk)
mean(chunk_na$tab$nas)

# check components of scale_ndvi

f1 <- function(x) floor(x/10)/1000
e <- t(apply(chunk,1,f1))
# everything will be converted into NA
length(which(e <= 0))
# there are no positive values at all in the floored data
length(which(e > 0))

f2 <- function(x) x - floor(x/10)*10 +1
r <- t(apply(chunk,1,f2))
# there some non-zero/non-NA values of r
mean(r,na.rm = TRUE)
# but none are as large as 4
r[which(r >= 4)]


# because all scaled values are <= 0, this operation 
# turns everything into NA
e = ifelse(e <= 0,NA,e)
e[which(r == 4 | r == 6)] = 0
e[which(r == 7)] = NA

# this logic returns ALL NA values 
e_na <- DescTools::CountCompCases(e)
mean(e_na$tab$nas)

# re-run to see where NAs come from
e <- t(apply(chunk,1,f1))
# none at this stage
e_na <- DescTools::CountCompCases(e)
mean(e_na$tab$nas)

# this turns everything into NA
# because there are no scaled values greater than 0
e[which(e > 0)]
e = ifelse(e <= 0,NA,e)
e_na <- DescTools::CountCompCases(e)
mean(e_na$tab$nas)

# this doesn't turn anything back to 0
# because all values of r < 4
e[which(r == 4 | r == 6)] = 0
e[which(r == 7)] = NA
e_na <- DescTools::CountCompCases(e)
mean(e_na$tab$nas)
