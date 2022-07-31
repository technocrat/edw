m <- matrix(NA,nrow=12,ncol=12)
colnames(m) <- colnames(y)
v1 <- 1:12
v2 <- 13:24
v3 <- 25:36
v4 <- 37:48
v5 <- 49:60
v6 <- 61:72
v7 <- 73:84
v8 <- 85:96
v9 <- 97:108
v10 <- 109:120
v11 <- 121:132
v12 <- 133:144
v <- list(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12)
for(i in seq_along(v)) m[i,] = v[[i]]
