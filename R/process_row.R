m <- matrix(NA,nrow=dim(leftover)[1],ncol=12)
# M <- Matrix(NA, dim(leftover)[1], 12) about the same time to produce nothing
colnames(m) <- colnames(receiver)
process_row <- function(x) {
  for(i in nrow(x)){
    # accumulate results in a vector
    v = vector(mode = "list", length = 12)
    v[1] = x[i,1] # pix value for this chunk
    # discard pix
    the_chunk = x[1,-1]
    the_chunk = scale_ndvi(the_chunk)
    v[11] = length(which(is.na(the_chunk)))
    v[12] = length(
      which(the_chunk - 
              floor(the_chunk/10)*10 + 1 == 4 | 
              the_chunk - 
              floor(the_chunk/10)*10 + 1 == 6))
    # 2022-07-30 14:18 -07:00 approx requires at least two NA
    if(v[12] > 1) the_chunk = interpolate(the_chunk)
    the_chunk = make_splines(the_chunk)
    v[2:8] = make_smooth(the_chunk)
    v[9] = make_sss(the_chunk)
    v[10] = L.JDAY
    v = unlist(t(v))
    m[i,] = unlist(t(v))
    return(m[i,])
  }
  rbind(m,m[i,])
  return(m)
}
# 0.6 sec
tic()
output <- process_row(leftover)
toc()
