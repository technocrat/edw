# scaling_test.R
# test effect of scaling algorithm
# on first 2000 pixel rows
# author: Richard Careaga
# Date: 2022-08-10

# LIBRARIES
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

mean(d,na.rm = TRUE)
min(d)
max(d)

# unscaled NDVI ranges from -10000 to +10000, so the most
# straighforward scaling is simply to multiply by 1e-4

simple <- d * 1e-4 
mean(simple,na.rm = TRUE)
min(simple)
max(simple)

# original dopixel code
# data.e<-floor(data[p,]/10)/1000
# my equivalent 
# e = floor(x/10)/1000

# very slightly less than simple
o <- floor(d/10)/1000
mean(o,na.rm = TRUE)
min(o,na.rm = TRUE)
max(o,na.rm = TRUE)

# original dopixel code
# data.e[data.e<=0]<-NA

o2 <- o
o2[o2 <= 0] <- NA
mean(o2,na.rm = TRUE)

# my equivalent
e = ifelse(o <= 0,NA,o)
min(e,na.rm = TRUE)
max(e,na.rm = TRUE)
mean(e, na.rm = TRUE)

# data.r
# original dopixel code
# data.r<-data[p,]-floor(data[p,]/10)*10 +1
# my equivalent
r <- d - floor(d/10)*10 + 1
min(r,na.rm = TRUE)
max(r,na.rm = TRUE)
mean(r, na.rm = TRUE)



# apply data.r to data.e
# original dopixel code
# data.e[which(data.r==4 | data.r==6)]<-0
# my equivalent
e[which(r == 4 | r == 6)] <- 0
min(e,na.rm = TRUE)
max(e,na.rm = TRUE)
mean(e, na.rm = TRUE)

# original dopixel code
# data.e[data.r==7]<-NA
e[which(r == 7)] <- NA
min(e,na.rm = TRUE)
max(e,na.rm = TRUE)
mean(e, na.rm = TRUE)
rowMeans(e,na.rm = TRUE) |> quantile(x = _, seq(0,1,0.1))

# function identical
E <- scale_ndvi(d)
rowMeans(E,na.rm = TRUE) |> quantile(x = _, seq(0,1,0.1))
# quantiles
qs <- t(apply(e,1,f))
lq <- qs[,1]
uq <- qs[,2]
am <- lq-uq