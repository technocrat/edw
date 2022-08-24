# determine how scale_ndvi should be applied
# to d or to d rowwise
x <- d
e = floor(x/10)/1000
r = x - floor(x/10)*10 +1
e = ifelse(e <= 0,NA,e)
e[which(r == 4 | r == 6)] = 0
e[which(r == 7)] = NA
mean(scale_ndvi(d),na.rm = TRUE)
mean(d, na.rm = TRUE)
mean(d,1,scale_ndvi,na.rm = TRUE)
mean(scale_ndvi(d[1,]), na.rm = TRUE)
mean(apply(d,1,scale_ndvi),na.rm = TRUE)
# produces inverted return matrix
g <- apply(d,1,scale_ndvi)
# this looks right
g <- t(apply(d,1,scale_ndvi))
rm(e,f,g,r,x)
