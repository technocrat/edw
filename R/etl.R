# etl.R
# extract, transform load csv file
# author: Richard Careaga
# Date: 2022-08-03

# populates a matrix, M, in sourcing script with an object,
# chunk, put into the environment by the sourcing
# script

M[,11] = length(which(is.na(chunk))) # NA.length
M[,12] = length(                     # ICE.length
    which(chunk - 
            floor(chunk/10)*10 + 1 == 4 | 
            chunk - 
            floor(chunk/10)*10 + 1 == 6))
# interpolate all rows
chunk = t(apply(chunk,1,interpolate))
begin = chunk[,1]
finish = chunk[,780]
chunk <- make_splines(chunk)
# chunk is used for find_mean_abs and make_sss
# chunk_ss is used for the smooth.spline objects
# filter is circular so take first value from data
chunk_ss <- chunk
chunk_ss[,1] = begin
# filter is circular so take last value from data
chunk_ss[,780] = finish
# TODO rep?
ss = smooth.spline(rep(JDAY.x,nrow(chunk)), t(chunk_ss), all.knots = TRUE)$y
# these need to be rowwise
M[,2] = quantile(ss, 0.1)             # lq
M[,3] = quantile(ss, 0.9)             # uq
M[,4] = M[,3] - M[,2]                 # amplitude
M[,5] = mean(ss)                      # mean.evi
M[,6] = sd(ss)                        # sd.evi
M[,7] = sum(ss)                       # sum.evi
M[,8] = apply(chunk,1,find_mean_abs)  # devi.ss
M[,9] = apply(chunk,1,make_sss)       # devi.sss
M[,10] = L.JDAY                       # L.JDAY


