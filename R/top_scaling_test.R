# scaling_test.R
# test effect of scaling algorithm
# author: Richard Careaga
# Date: 2022-08-10

# LIBRARES
library(here)

# DATA
# `d` is equivalent to `data` in original scripts
# data, df, c, t and other built-in names should
# be avoided because some operations give priority
# in namespace to the built-in
d <- read.csv(here("obj/chunks/chunk_1.csv"))

# convert to matrix

d <- as.matrix(d)

dim(d)

# remove 'pix` identifier column
d <- d[,-1]
dim(d)

# average value across all dates and rows
# shows that there are no NA values in data

mean(d)
min(d)
max(d)

# unscaled NDVI ranges from -10000 to +10000, so the most
# straighforward scaling is simply to multiply by 1e-4

simple <- d * 1e-4 
mean(simple)
min(simple)
max(simple)

# original dopixel code
# data.e<-floor(data[p,]/10)/1000
# my equivalent 
# e = floor(x/10)/1000

o <- floor(d/10)/1000
mean(o)
min(o)
max(o)

# rowwise produces same result
f <- function(x) floor(x/10)/1000
orw <- apply(d,1,f)
mean(orw)
min(orw)
max(orw)
# original dopixel code
# data.e[data.e<=0]<-NA

o2 <- o
o2[o2 <= 0] <- NA
mean(o2,na.rm = TRUE)

# my equivalent
e = ifelse(o <= 0,NA,o)
min(e)
max(e)
mean(e, na.rm = TRUE)

# data.r
# original dopixel code
# data.r<-data[p,]-floor(data[p,]/10)*10 +1
# my equivalent
r <- d - floor(d/10)*10 + 1
min(r)
max(r)
mean(r, na.rm = TRUE)

# rowwise produces same result
f <- function(x) d[x,] - floor(d[x,]/10)*10 + 1
rrw <- apply(d,1,f)
rrw <- d[1:2000,] - floor(d[1:2000,]/10)*10 + 1
min(rrw)
max(rrw)
mean(rrw, na.rm = TRUE)

# apply data.r to data.e
# original dopixel code
# data.e[which(data.r==4 | data.r==6)]<-0
# my equivalent
e[which(r == 4 | r == 6)] <- 0
min(e)
max(e)
mean(e, na.rm = TRUE)

# original dopixel code
# data.e[data.r==7]<-NA
# won't do anything because r has no values of 7
e[which(r == 7)] <- NA
min(e)
max(e)
mean(e, na.rm = TRUE)
rowMeans(e,na.rm = TRUE) |> quantile(x = _, seq(0,1,0.1))

