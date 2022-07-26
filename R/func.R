# func.R
# functions used in scripts under this directory
# author: Richard Careaga
# Date: 2022-07-20

# how many NA data points on scaled data 

count_na <- function(x,y = chunk_size) {
  v <- vector(length = chunk_size)
  for(i in 1:y) v[i] = length(which(is.na(x)))
  return(v)
}


# takes a row vector of the scaled data and returns a scalar

find_frozen <- function(x,y = chunk_size){
  # how many possible ICE/SNOW data points
  v <- vector()
  # remove pix column
  x = x[,-1]
  for(i in 1:y) {
    #q = x[i,] - floor(unlist(x[i,])/10)*10 + 1
    v[i] = length(which((unlist(x[i,]) - 
                           floor(unlist(x[i,])/10)*10 + 1) == 4 | 
                        (unlist(x[i,]) - 
                           floor(unlist(x[i,])/10)*10 + 1) == 6))  
  }
  return(v)
}

# interpolate missing values of scaled data
 
interpolate <- function(x,y=chunk_size) {
  m = matrix(0,chunk_size,col_width)
  for(i in 1:chunk_size) {
    m[i,] = approx(JDAY.x,unlist(x[i,-1]), xout = JDAY.x,  rule = 2)$y
  }
  return(m)
}

# takes scaled data

make_splines <- function(x) {
  # set aside record id
  pix = x[,1]
  x   = x[,-c(1)]
  # get a filtered running mean:
  splines = filter(x,
        filter = rep(1 / 3,3),
        method = "convolution",
        sides = 2, # 2022-07-21 00:20 -07:00 "side" in the original
        circular = TRUE
    )
  # reunite with pix
  splines = cbind(pix,splines)
  return(splines)
} 


make_smooth <- function(x){
  # adjust col width to reflect removal of first and last values
  #ss = matrix(0,chunk_size,col_width - 2) 
  for(i in 1:chunk_size)  {
    # TODO fix hardwired 781
    ss = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                       x[i,-c(1:2,781)])$y
    # TODO: need to capture pix, not i
    pix = i
    lq = quantile(ss, 0.1)
    uq = quantile(ss, 0.9)
    amplitude = uq - lq
    # TODO: write to receiver object
    # the mean evi of the time series
    mean.evi = mean(ss)
    # the standard deviation of the time series
    sd.evi = sd(ss)
    # sum evi
    sum.evi = sum(ss)
    # deviation between less smooth spline and data
    devi.ss <- mean(unlist(abs(s[i,-c(1:2,781)]-ss))) 
    results = c(pix,lq,uq,amplitude,mean.evi,sd.evi,sum.evi,devi.ss)
    names(results) <- c("pix","lq","uq","amplitude","mean.evi","sd.evi","sum.evi","devi.ss")
    return(results)
  }
}

# takes the scaled data

make_sss <- function(x){
  # adjust col width to reflect removal of first and last values
  #sss <- matrix(0,chunk_size,col_width - 2) 
  # return object
  v = vector()
  for(i in 1:chunk_size)  {
    # TODO fix hardwired 781
    sss_both = smooth.spline(JDAY.x[-c(2,col_width - 1)],
                         x[i,-c(1:2,781)],
                         spar = 0.3)
    sss.x = sss_both$x
    sss.y = sss_both$y
    d1 = predict(sss_both, sss.x, deriv = 1)$y
    d1 = (d1 - mean(d1)) / sd(d1)
    d2 = predict(sss_both, sss.x, deriv = 2)$y
    d2 = (d2 - mean(d2)) / sd(d2)
    # deviation between smoother spline and data
    devi.sss = mean(unlist(abs(x[i,-c(1:2,781)] - sss.y)))  
    return(list(d1,d2,devi.sss))
    }
}

# save the current chunk's pix column

save_pix <- function(x) x[,1]

# scale source data

scale_ndvi <- function(x) {
  # save out record id, pix
  pix = x[,1]
  x   = x[,-c(1)]
  e = floor(x/10)/1000
  r = x - floor(x/10)*10 +1
  ifelse(e <= 0,NA,0)
  e[which(r == 4 | r == 6)] = 0
  e[which(r == 7)] = NA
  # restore pix
  e = cbind(pix,e)
  return(e)
}

