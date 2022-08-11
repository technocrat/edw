# do this in /home, on /media IO operations slow everything down
# presumable from reading and writing data between different disk partitions


library(raster)
library(ncdf4)
library(ncdf4.helpers)
library(bigmemory)
library(foreach)
library(doMC)
library(doParallel)
#library(conflicted)

#conflict_prefer("loadings", "stats")
cores<-30
registerDoMC(cores=cores)
#cores<-getDoParWorkers()


#--------------------------------------------------------------------------------------------------
rad <- function(degree) (degree * pi)/180

deg <- function(radian) (radian * 180)/pi

mycv <- function( x ) (var(x)^0.5) / mean(x)

circ.mean <- function(x) {
  x<-rad(x*(360/365))
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  circmean <- deg(circmean)
  if(circmean<0) circmean<-365+circmean
  circmean}
  
trap1<-function( x, a,  b) #trapezoid1
{
  pmax( pmin( (x-a)/(b-a), 1 ), 0 )
}
  
setwd("~/fbiome_data")
#--------------------------------------------------------------------------------------------------
# set day counters
data<-attach.big.matrix("ndvi.desc")
tail(data)
date<-as.POSIXlt(substr(colnames(data),2,11), tz="GMT", format="%Y_%m_%d")
YEAR<-date$year #year
JDAY<-date$yday #julian day
# Ugly way to set leapyear dates to same Julian day as non-leapyear dates (for peak.day in dopixel)
JDAY<-ifelse(c(JDAY==67|JDAY==82|JDAY==98|JDAY==113|JDAY==128|JDAY==143|JDAY==159|JDAY==174|JDAY==189|
                 JDAY==204|JDAY==220|JDAY==235|JDAY==251|JDAY==266|JDAY==281|JDAY==296|JDAY==312|JDAY==327|
                 JDAY==342|JDAY==357),JDAY-1,JDAY)
JDAY.x<-JDAY+(YEAR-min(YEAR))*365



#---------------------Merra-2 clear sky netcdf
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(chron)

library(bigmemory)
library(biganalytics)
library(foreach)
library(doMC)



#######MERRA DATA IN


etemp<-raster("MERRA2_400.tavgM_2d_rad_Nx.1980_2015.nc4",varname="SWGNTCLRCLN") #a single years data
ncin <- nc_open("MERRA2_400.tavgM_2d_rad_Nx.1980_2015.nc4")
print(ncin)

names(ncin$var)


lon <- ncvar_get(ncin, "lon")
nlon <- dim(lon)
head(lon)
lat <- ncvar_get(ncin, "lat", verbose = F)
nlat <- dim(lat)
head(lat)


# Get the time variable and its attributes using the ncvar_get() and 
# ncatt_get() functions, and also get the number of times using the dim() function.
print(c(nlon, nlat))
t <- ncvar_get(ncin, "time")
tunits <- ncatt_get(ncin, "time", "units")
nt <- dim(t)

head(dates(t))

mydates = t
mydates = as.POSIXct(t*60,origin = "1980-01-01",tz="GMT")
head(mydates)

ct<-mydates

#Get the variable and its attributes, and verify the size of the array.
#"SWGDN"

dname<-"SWGNTCLRCLN" #"SWGDN"
#"SWGNTCLRCLN"
rad.array <- ncvar_get(ncin, dname)
dlname <- ncatt_get(ncin, dname, "long_name")
dunits <- ncatt_get(ncin, dname, "units")
fillvalue <- ncatt_get(ncin, dname, "_FillValue")
dim(rad.array)

rad.array<-rad.array[,,1:420] #take last year out of array because of some time issues with time, maybe data is incomplete anyhow
ct<-ct[1:420]


#ncvar_get
#Close the NetCDF file using the nc_close() function.
nc_close(ncin)

m.date<-as.POSIXlt(ct, tz="GMT", format="%Y.%m.%d")
m.YEAR<-m.date$year #year as.integer(format(m.date,"%Y")) #year
m.JDAY<-m.date$yday +15  #julian day #push to midlle of month and not 1 of month to better represent montly mean
# Ugly way to set leapyear dates to same Julian day as non-leapyear dates (for peak.day in dopixel)
m.JDAY<-ifelse(c(m.JDAY==67|m.JDAY==82|m.JDAY==98|m.JDAY==113|m.JDAY==128|m.JDAY==143|m.JDAY==159|m.JDAY==174|m.JDAY==189|
                 m.JDAY==204|m.JDAY==220|m.JDAY==235|m.JDAY==251|m.JDAY==266|m.JDAY==281|m.JDAY==296|m.JDAY==312|m.JDAY==327|
                 m.JDAY==342|m.JDAY==357),m.JDAY-1,m.JDAY)
m.JDAY.x<-m.JDAY+(m.YEAR-81)*365 #NOTE use YEAR=1981 for min here not m.YEAR!


#xy are the xy coords of the NDVI data
xy<-attach.big.matrix("xy.desc")
xy<-SpatialPoints(cbind(xy[,1],xy[,2]),proj4string=crs(etemp))


#print("#make a raster with sol's grid specifications (r.index)")
# etemp<-raster("MERRA2_400.tavgM_2d_rad_Nx.1980_2015.nc4",varname="SWGNTCLRCLN") #a single years data
# r.index<-raster(etemp[[1]])
# #fill this raster with integers 1:ncell(r.index)
# r.index<-setValues(r.index,1:ncell(r.index))
# #extract cell numbers e.g. sol.row<-extract(r.index,coordinates(xy))
# txy<-cbind(coordinates(xy)[,1],coordinates(xy)[,2])
# sol.row<-extract(r.index,txy)
# 

txy<-cbind(coordinates(xy)[,1],coordinates(xy)[,2])
txy[,2]<-round(coordinates(xy)[,2]*(2/1)) / (2/1) #lat y row 361
txy[,1]<-round(coordinates(xy)[,1]*(8/5)) / (8/5) #lon x col 576
txy[,1][txy[,1]==180]<- -180 #no 180 only -180

# #testing coord lookup system
# p<-2051276
# p<-1895392
# p<-2074747
# p<-725893
# p<-1648046
# dfss<-data.frame(d=m.JDAY.x,sr=rad.array[ which(lon==txy[p,1]),which(lat==txy[p,2]),]) #dim(rad.array) 576 361 420
# plot(dfss$d[1:24],dfss$sr[1:24],type="l");dev.off()
# coordinates(xy)[p,]
# txy[p,]



#######2018 MERRA DATA IN
rm(ncin)


ncin.land <- nc_open("MERRA2_400.tavgM_2d_lnd_Nx.1980_2018.nc4")

print(ncin.land)

names(ncin.land$var)
# "EVLAND"      "GWETPROF"    "GWETROOT"    "PRECTOTLAND" "RZMC"       
# "SNODP"       "SWLAND"      "TSOIL1"  

etemp.land<-raster("tmp.MERRA2_400.tavgM_2d_lnd_Nx.201507.nc4.nc",varname="SWLAND") #a single years data


land.lon <- ncvar_get(ncin.land, "lon")
land.nlon <- dim(land.lon)
head(land.lon)
land.lat <- ncvar_get(ncin.land, "lat", verbose = F)
land.nlat <- dim(land.lat)
head(land.lat)


# Get the time variable and its attributes using the ncvar_get() and 
# ncatt_get() functions, and also get the number of times using the dim() function.
print(c(land.nlon, land.nlat))
t <- ncvar_get(ncin.land, "time")
tunits <- ncatt_get(ncin.land, "time", "units")
nt <- dim(t)

head(dates(t))

mydates = t
mydates = as.POSIXct(t*60,origin = "1980-01-01",tz="GMT")
head(mydates)

ct<-mydates

#Get the variable and its attributes, and verify the size of the array.
#"SWGDN"

extract.var<-function(dname) {
#dname<-"SWLAND" 
the.array <- ncvar_get(ncin.land, dname)
dlname <- ncatt_get(ncin.land, dname, "long_name")
dunits <- ncatt_get(ncin.land, dname, "units")
fillvalue <- ncatt_get(ncin.land, dname, "_FillValue")
the.array[the.array == fillvalue$value] <- NA
#image(the.array[,,100])

dim(the.array)

the.array<-the.array[,,1:420] #take last year out of array because of some time issues with time, maybe data is incomplete anyhow
return(the.array)
}

SWLAND.array<-extract.var("SWLAND")
RZMC.array<-extract.var("RZMC")



ct<-ct[1:420]


#ncvar_get
#Close the NetCDF file using the nc_close() function.
nc_close(ncin.land)

merra.date<-as.POSIXlt(ct, tz="GMT", format="%Y.%merra.%d")
merra.YEAR<-merra.date$year #year as.integer(format(merra.date,"%Y")) #year
merra.JDAY<-merra.date$yday +15  #julian day #push to midlle of month and not 1 of month to better represent montly mean
# Ugly way to set leapyear dates to same Julian day as non-leapyear dates (for peak.day in dopixel)
merra.JDAY<-ifelse(c(merra.JDAY==67|merra.JDAY==82|merra.JDAY==98|merra.JDAY==113|merra.JDAY==128|merra.JDAY==143|merra.JDAY==159|merra.JDAY==174|merra.JDAY==189|
                 merra.JDAY==204|merra.JDAY==220|merra.JDAY==235|merra.JDAY==251|merra.JDAY==266|merra.JDAY==281|merra.JDAY==296|merra.JDAY==312|merra.JDAY==327|
                 merra.JDAY==342|merra.JDAY==357),merra.JDAY-1,merra.JDAY)
merra.JDAY.x<-merra.JDAY+(merra.YEAR-81)*365 #NOTE use YEAR=1981 for min here not merra.YEAR!

 

 #xy are the xy coords of the NDVI data
merra.xy<-attach.big.matrix("xy.desc")
merra.xy<-SpatialPoints(cbind(merra.xy[,1],merra.xy[,2]),proj4string=crs(etemp.land))

merra.txy<-cbind(coordinates(merra.xy)[,1],coordinates(merra.xy)[,2])
merra.txy[,2]<-round(coordinates(merra.xy)[,2]*(2/1)) / (2/1) #lat y row 361
merra.txy[,1]<-round(coordinates(merra.xy)[,1]*(8/5)) / (8/5) #lon x col 576
merra.txy[,1][txy[,1]==180]<- -180 #no 180 only -180




####### TEMPERATURE DATA IN
ncin <- nc_open("cru_ts3.22.1971.2013.tmp.dat.nc")
print(ncin)
tmp.lon <- ncvar_get(ncin, "lon")
nlon <- dim(tmp.lon)
head(tmp.lon)
tmp.lat <- ncvar_get(ncin, "lat", verbose = F)
nlat <- dim(tmp.lat)
head(tmp.lat)

# Get the time variable and its attributes using the ncvar_get() and 
# ncatt_get() functions, and also get the number of times using the dim() function.
print(c(nlon, nlat))
tmp.t <- ncvar_get(ncin, "time")
tmp.tunits <- ncatt_get(ncin, "time", "units")
tmp.nt <- dim(tmp.t)

tmp.ct<-as.Date(tmp.t,origin="1900-01-01") 

#Get the variable and its attributes, and verify the size of the array.
dname<-"tmp"
tmp.array <- ncvar_get(ncin, dname)
dlname <- ncatt_get(ncin, dname, "long_name")
dunits <- ncatt_get(ncin, dname, "units")
fillvalue <- ncatt_get(ncin, dname, "_FillValue")
dim(tmp.array)
#ncvar_get
#Close the NetCDF file using the nc_close() function.
nc_close(ncin)

tmp.m.date<-as.POSIXlt(tmp.ct, tz="GMT", format="%Y.%m.%d")
tmp.m.YEAR<-tmp.m.date$year #year as.integer(format(m.date,"%Y")) #year
tmp.m.JDAY<-tmp.m.date$yday +15  #julian day #push to midlle of month and not 1 of month to better represent montly mean
# Ugly way to set leapyear dates to same Julian day as non-leapyear dates (for peak.day in dopixel)
tmp.m.JDAY<-ifelse(c(tmp.m.JDAY==67|tmp.m.JDAY==82|tmp.m.JDAY==98|tmp.m.JDAY==113|tmp.m.JDAY==128|tmp.m.JDAY==143|tmp.m.JDAY==159|tmp.m.JDAY==174|tmp.m.JDAY==189|
                 tmp.m.JDAY==204|tmp.m.JDAY==220|tmp.m.JDAY==235|tmp.m.JDAY==251|tmp.m.JDAY==266|tmp.m.JDAY==281|tmp.m.JDAY==296|tmp.m.JDAY==312|tmp.m.JDAY==327|
                 tmp.m.JDAY==342|tmp.m.JDAY==357),tmp.m.JDAY-1,tmp.m.JDAY)
tmp.m.JDAY.x<-tmp.m.JDAY+(tmp.m.YEAR-81)*365 #NOTE use YEAR=1981 for min here not m.YEAR!


#this already done above
# xy<-attach.big.matrix("xy.desc")
# xy<-SpatialPoints(cbind(xy[,1],xy[,2]),proj4string=crs(etemp))


# tmp.txy<-cbind(coordinates(xy)[,1],coordinates(xy)[,2])
# tmp.txy[,2]<-round(coordinates(xy)[,2]*(2/1)) / (2/1) + 0.25 #lat y row 360
# tmp.txy[,1]<-round(coordinates(xy)[,1]*(2/1)) / (2/1) + 0.25 #lon x col 720
# tmp.txy[,1][txy[,1]>180]<- -179.25 #no 180.75 only -179.75

 
 xy2<-attach.big.matrix("xy.desc")
 xr<-floor(xy2[,1]*2)/2+0.25
 yr<-floor(xy2[,2]*2)/2+0.25
 tmp.txy<-cbind(xr,yr)


## SOIL MOISTURE IN


ncin <- nc_open("soilw.mon.mean.v2.nc")

print(ncin)
sm.lon <- ncvar_get(ncin, "lon")
sm.nlon <- dim(sm.lon)
sm.lon[sm.lon>180]<-sm.lon[sm.lon>180] - 360
summary(sm.lon)

sm.lat <- ncvar_get(ncin, "lat", verbose = F)
sm.nlat <- dim(sm.lat)
#sm.lat <- -1*sm.lat #seems wrong way round! as revealed by image(soilw.array[,,1])
head(sm.lat)



# Get the time variable and its attributes using the ncvar_get() and 
# ncatt_get() functions, and also get the number of times using the dim() function.
print(c(sm.nlon, sm.nlat))
sm.t <- ncvar_get(ncin, "time")
sm.tunits <- ncatt_get(ncin, "time", "units")
sm.nt <- dim(sm.t)

sm.ct<-as.Date(sm.t,origin="1800-01-01") 

#Get the variable and its attributes, and verify the size of the array.
dname<-"soilw"
soilw.array <- ncvar_get(ncin, dname)
dlname <- ncatt_get(ncin, dname, "long_name")
dunits <- ncatt_get(ncin, dname, "units")
fillvalue <- ncatt_get(ncin, dname, "_FillValue")
dim(soilw.array)

#Get the global attributes.
title <- ncatt_get(ncin, 0, "title")
institution <- ncatt_get(ncin, 0, "institution")
datasource <- ncatt_get(ncin, 0, "source")
references <- ncatt_get(ncin, 0, "references")
history <- ncatt_get(ncin, 0, "history")
Conventions <- ncatt_get(ncin, 0, "Conventions")

#ncvar_get
#Close the NetCDF file using the nc_close() function.
nc_close(ncin)

# Replace NetCDF fillvalues with R NAs
soilw.array[soilw.array == fillvalue$value] <- NA


 soilw.array2<-soilw.array[,nlat:1,1]
 image(soilw.array2[,])
 image(soilw.array[,nlat:1,1])
 dev.off()
 
 soilw.array<-soilw.array[,nlat:1,]
class(soilw.array)

sm.m.date<-as.POSIXlt(sm.ct, tz="GMT", format="%Y.%m.%d")
sm.m.YEAR<-sm.m.date$year #year as.integer(format(m.date,"%Y")) #year
sm.m.JDAY<-sm.m.date$yday +15  #julian day #push to midlle of month and not 1 of month to better represent montly mean
# Ugly way to set leapyear dates to same Julian day as non-leapyear dates (for peak.day in dopixel)
sm.m.JDAY<-ifelse(c(sm.m.JDAY==67|sm.m.JDAY==82|sm.m.JDAY==98|sm.m.JDAY==113|sm.m.JDAY==128|sm.m.JDAY==143|sm.m.JDAY==159|sm.m.JDAY==174|sm.m.JDAY==189|
                 sm.m.JDAY==204|sm.m.JDAY==220|sm.m.JDAY==235|sm.m.JDAY==251|sm.m.JDAY==266|sm.m.JDAY==281|sm.m.JDAY==296|sm.m.JDAY==312|sm.m.JDAY==327|
                 sm.m.JDAY==342|sm.m.JDAY==357),sm.m.JDAY-1,sm.m.JDAY)
sm.m.JDAY.x<-sm.m.JDAY+(sm.m.YEAR-81)*365 #NOTE use YEAR=1981 for min here not m.YEAR!


#this already done above
# xy<-attach.big.matrix("xy.desc")
# xy<-SpatialPoints(cbind(xy[,1],xy[,2]),proj4string=crs(etemp))


# sm.txy<-cbind(coordinates(xy)[,1],coordinates(xy)[,2])
# sm.txy[,2]<-round(coordinates(xy)[,2]*(2/1)) / (2/1) + 0.25 #lat y row 360
# sm.txy[,1]<-round(coordinates(xy)[,1]*(2/1)) / (2/1) + 0.25 #lon x col 720
# sm.txy[,1][txy[,1]>180]<- -179.25 #no 180.75 only -179.75

 xy2<-attach.big.matrix("xy.desc")
 xr<-floor(xy2[,1]*2)/2+0.25
 yr<-floor(xy2[,2]*2)/2+0.25
# xr <- 180-xr #convert from netcdf 0 to 360
 sm.txy<-cbind(xr,yr)

#setwd("/home/higgins")
#setwd("/home/higgins/BiomeChange/working")



#-END STEVE-SETUP SOLRAD----------------------------------------------------------------------------------------------


p<-270000


rm(date)


#--------------------------------------------------------------------------------------------------
# create big matrix for stats
# colnames for bigmatrix
names<-c("lq", "uq", "mean.evi", "sd.evi", "sum.evi", "amplitude",
         "peak.day", "trough.day",
         "devi.ss", "devi.sss", "L.JDAY",
         "NA.length", "ICE.length",
         "cor.photofac.all","cor.moistfac.all","cor.radfac.all", "cor.evi.all","cv.gpp.all","cv.evi.all",
         "cor.vpi.photofac.all","cor.vpi.moistfac.all","cor.vpi.radfac.all", "cor.vpi.evi.all","cv.vpi.all",
         paste("td", 1:33, sep="."), paste("td.x", 1:33, sep="."), paste("td.evi", 1:33, sep="."),
         paste("pd", 1:32, sep="."), paste("pd.x", 1:32, sep="."), paste("pd.evi", 1:32, sep="."),
         paste("elon.m", 1:32, sep="."), paste("elon.m.x", 1:32, sep="."),
         paste("eoff.m", 1:32, sep="."), paste("eoff.m.x", 1:32, sep="."),
         paste("elon.f", 1:32, sep="."), paste("elon.f.x", 1:32, sep="."),
         paste("eoff.f", 1:32, sep="."), paste("eoff.f.x", 1:32, sep="."),
         paste("elon.l", 1:32, sep="."), paste("elon.l.x", 1:32, sep="."),
         paste("eoff.l", 1:32, sep="."), paste("eoff.l.x", 1:32, sep="."),
         paste("sum.evi.yr", 1:32, sep="."), paste("amp", 1:32, sep="."),
         paste("gsl", 1:32, sep="."), paste("gsl.peak", 1:32, sep="."), paste("gsl.long", 1:32, sep="."),
         paste("elon.f.evi", 1:32, sep="."), paste("elon.m.evi", 1:32, sep="."), paste("elon.l.evi", 1:32, sep="."),
         paste("eoff.f.evi", 1:32, sep="."), paste("eoff.m.evi", 1:32, sep="."), paste("eoff.l.evi", 1:32, sep="."),
         paste("sum.gpp.yr", 1:32, sep="."), #paste("td.gpp.yr", 1:32, sep="."),paste("pd.gpp.yr", 1:32, sep="."),
         paste("sum.vpi.yr", 1:32, sep="."),
         paste("sum.sol.yr", 1:32, sep="."),
         paste("td.gpp", 1:32, sep="."),paste("pd.gpp", 1:32, sep="."),
         paste("min.gpp", 1:32, sep="."),paste("max.gpp", 1:32, sep="."),
         paste("amp.gpp", 1:32, sep="."),
         paste("tdg", 1:32, sep="."),paste("pdg", 1:32, sep="."),
         paste("tdg.x", 1:32, sep="."),paste("pdg.x", 1:32, sep="."),
         paste("cor.photofac",1:32, sep="."),paste("cor.moistfac",1:32, sep="."),paste("cor.radfac",1:32, sep="."),
         paste("cor.evi",1:32, sep="."), paste("cv.gpp",1:32, sep="."), paste("cv.evi",1:32, sep="."), #2018
         paste("cor.vpi.photofac",1:32, sep="."),paste("cor.vpi.moistfac",1:32, sep="."),paste("cor.vpi.radfac",1:32, sep="."),
         paste("cor.vpi.evi",1:32, sep="."), paste("cv.vpi",1:32, sep=".") #2018
         ) 

#pd.gpp and pd.gpp.yr are duplicates
#td.gpp and td.gpp.yr are duplicates
         
stats<-big.matrix(nrow=nrow(data), ncol=length(names), type="double",
                  dimnames=list(row.names=NULL, col.names=names),
                  backingfile="stats.bin", descriptorfile="stats.desc")

                  
#Use shared descriptor, we can load the data from disk with
#attach.big.matrix function.
                  
#--------------------------------------------------------------------------------------------------
print("# run dopixel for each row in the big matrix")
source("dopixel_3g_v2_01.R")
library("pls")


# cores<-32
# registerDoMC(cores=cores)


#source("dopixel_3g_v2_01.R")
statsbio<-attach.big.matrix("stats.desc")  
registerDoParallel(30)
tot.count<-120
start<-round( seq(1, nrow(data), length.out=tot.count+1), digits=0)
end<-c(start[(2):tot.count]-1,start[tot.count+1])
for(count in 1:120) {
gc()
print(count)
foreach(p = start[count]:end[count]) %dopar% {
    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")
  }
}  
stopImplicitCluster()

sessionInfo()
str()

saveRDS(merra.txy, "merra.txy.rds")
saveRDS(sm.txy, "sm.txy.rds")
saveRDS(soilw.array2, "soilw.array2.rds")
saveRDS(tmp.txy, "tmp.txy.rds")
saveRDS(txy, "txy.rds")

saveRDS(Conventions, "Conventions.rds")
saveRDS(ct, "ct.rds")
saveRDS(data, "data.rds")
saveRDS(datasource, "datasource.rds")
saveRDS(dlname, "dlname.rds")
saveRDS(dunits, "dunits.rds")
saveRDS(etemp, "etemp.rds")
saveRDS(etemp.land, "etemp.land.rds")
saveRDS(fillvalue, "fillvalue.rds")
saveRDS(history, "history.rds")
saveRDS(institution, "institution.rds")
saveRDS(JDAY, "JDAY.rds")
saveRDS(JDAY.x, "JDAY.x.rds")
saveRDS(land.lat, "land.lat.rds")
saveRDS(land.lon, "land.lon.rds")
saveRDS(land.nlat, "land.nlat.rds")
saveRDS(land.nlon, "land.nlon.rds")
saveRDS(lon, "lon.rds")
saveRDS(lat, "lat.rds")
saveRDS(m.date, "m.date.rds")
saveRDS(m.JDAY, "m.JDAY.rds")
saveRDS(m.JDAY.x, "m.JDAY.x.rds")
saveRDS(m.YEAR, "m.YEAR.rds")
saveRDS(merra.date, "merra.date.rds")
saveRDS(merra.JDAY, "merra.JDAY.rds")
saveRDS(merra.JDAY.x, "merra.JDAY.x.rds")
saveRDS(merra.xy, "merra.xy.rds")
saveRDS(merra.YEAR, "merra.YEAR.rds")
saveRDS(mydates, "mydates.rds")
saveRDS(names, "names.rds")
saveRDS(ncin, "ncin.rds")
saveRDS(ncin.land, "ncin.land.rds")
saveRDS(nlon, "nlon.rds")
saveRDS(nlat, "nlat.rds")
saveRDS(nt, "nt.rds")
saveRDS(rad.array, "rad.array.rds")
saveRDS(references, "references.rds")
saveRDS(RZMC.array, "RZMC.array.rds")
saveRDS(sm.ct, "sm.ct.rds")
saveRDS(sm.lon, "sm.lon.rds")
saveRDS(sm.lat, "sm.lat.rds")
saveRDS(sm.m.date, "sm.m.date.rds")
saveRDS(sm.m.JDAY, "sm.m.JDAY.rds")
saveRDS(sm.m.JDAY.x, "sm.m.JDAY.x.rds")
saveRDS(sm.m.YEAR, "sm.m.YEAR.rds")
saveRDS(sm.nlon, "sm.nlon.rds")
saveRDS(sm.nlat, "sm.nlat.rds")
saveRDS(sm.nt, "sm.nt.rds")
saveRDS(sm.t, "sm.t.rds")
saveRDS(sm.tunits, "sm.tunits.rds")
saveRDS(soilw.array, "soilw.array.rds")
saveRDS(SWLAND.array, "SWLAND.array.rds")
saveRDS(t, "t.rds")
saveRDS(title, "title.rds")
saveRDS(tmp.array, "tmp.array.rds")
saveRDS(tmp.ct, "tmp.ct.rds")
saveRDS(tmp.lon, "tmp.lon.rds")
saveRDS(tmp.lat, "tmp.lat.rds")
saveRDS(tmp.m.date, "tmp.m.date.rds")
saveRDS(tmp.m.JDAY, "tmp.m.JDAY.rds")
saveRDS(tmp.m.JDAY.x, "tmp.m.JDAY.x.rds")
saveRDS(tmp.m.YEAR, "tmp.m.YEAR.rds")
saveRDS(tmp.nt, "tmp.nt.rds")
saveRDS(tmp.t, "tmp.t.rds")
saveRDS(tmp.tunits, "tmp.tunits.rds")
saveRDS(tunits, "tunits.rds")
saveRDS(xr, "xr.rds")
saveRDS(xy, "xy.rds")
saveRDS(xy2, "xy2.rds")
saveRDS(YEAR, "YEAR.rds")
saveRDS(yr, "yr.rds")

##--1 Aug up to here
#c("bigmemory", "biganalytics", "biglm", "stats","utils", "methods", "pls")
#class(data)
registerDoParallel(30)

foreach(p = 1:nrow(data), .packages= "bigmemory") %dopar% {
    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")
    print(p)
  }
#first attempt got to at least 1 400 000 (p was 172 1411 on some workers)
stopImplicitCluster()




p<-start[12]+17520-1

for(p in start[count]:(start[count]+32*100))  {
    print(p)
    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")
  }


  
  
########################
sj<-seq(1,cores*10,by=cores)
ej<-seq(60,cores*10,by=cores)


for(count in 1:10) {  
print(c(count,sj[count],ej[count]))
  for(j in 1:5) print(c(j,start[j],end[j]))
}

foreach(j = sj[count]:ej[count]) %dopar% {
  print(c(j,start[j],end[j]))
#  for(p in start[j]:end[j]){
#    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")

#    print(p)
#  }
  }
}


for(count in 1:10) {  
print(count)
foreach(j = sj[count]:ej[count]) %dopar% {
  for(p in start[j]:end[j]){
    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")
#    print(p)
  }}
}
  

  

  
start<-round( seq(1, nrow(data), length.out=cores+1), digits=0)
end<-c(start[2:cores]-1,start[cores+1])

lon<-round(lon,3)
lat<-round(lat,3)
tmp.lon<-round(tmp.lon,3)
tmp.lat<-round(tmp.lat,3)

# 
# # For relatively small tasks it is more efficient to assign a for() loop to each worker
# # (letting foreach assign each task to a worker wastes significant overhead time)
# foreach(j = 1:cores) %dopar% {
#   for(p in start[j]:end[j]){
#     stats[p,] <- dopixel(p)
#   }}
#   

  
# 
start<-round( seq(1, nrow(data), length.out=cores*10+1), digits=0)
end<-c(start[2:(cores*10)]-1,start[cores*10+1])

sj<-seq(1,cores*10,by=cores)
ej<-seq(60,cores*10,by=cores)

for(count in 1:10) {  
print(count)
foreach(j = sj[count]:ej[count]) %dopar% {
  for(p in start[j]:end[j]){
    stats[p,] <- dopixel(p,PLOT=FALSE,REGION="NULL")
#    print(p)
  }}
}
# 
# start<-round( seq(nrow(data)/2+1, nrow(data), length.out=cores+1), digits=0)
# end<-c(start[2:cores]-1,start[cores+1])
# foreach(j = 1:cores) %dopar% {
#   for(p in start[j]:end[j]){
#     stats[p,] <- dopixel(p)
#   }}
#   
  
  


 
source("dopixel_3g_v2_01.R")

#stats[p,]<-stats.ret

 #amazon
 aw<-which ( coordinates(xy)[,2] < -5 & coordinates(xy)[,2] > -6 & coordinates(xy)[,1] < -62 & coordinates(xy)[,1] > -63)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Amazon_")
  }

 
 #pacific northwest 48, -120 forest near seattle
 aw<-which ( coordinates(xy)[,2] < 50 & coordinates(xy)[,2] > 49 & coordinates(xy)[,1] < -120 & coordinates(xy)[,1] > -121)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="PacificNW_")
  }
  
  
#florida/georgia  30.6, -82.4 
 aw<-which ( coordinates(xy)[,2] < 31 & coordinates(xy)[,2] > 29 & coordinates(xy)[,1] < -81 & coordinates(xy)[,1] > -83)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Florida_")
  }
    
  
 #canadian boreal 51 or 52.9, -94 or 96.7 
 aw<-which ( coordinates(xy)[,2] < 54 & coordinates(xy)[,2] > 52 & coordinates(xy)[,1] < -97 & coordinates(xy)[,1] > -98)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Boreal-Can_")
  }
  
 #russian boreal 60, 100.6
 aw<-which ( coordinates(xy)[,2] < 61 & coordinates(xy)[,2] > 59 & coordinates(xy)[,1] < 100 & coordinates(xy)[,1] > 99)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Boreal-Rus_")
  }
  
  

 
 #chiloe
 aw<-which ( coordinates(xy)[,2] < -41 & coordinates(xy)[,2] > -43 & coordinates(xy)[,1] < -73 & coordinates(xy)[,1] > -74)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Chiloe_")
  }
 
 
 #kruger
 aw<-which ( coordinates(xy)[,2] < -20 & coordinates(xy)[,2] > -21 & coordinates(xy)[,1] < 31 & coordinates(xy)[,1] > 30)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Kruger_")
  }
 
 #mali sahara 19 -3.5
 aw<-which ( coordinates(xy)[,2] < 20 & coordinates(xy)[,2] > 19 & coordinates(xy)[,1] < -3 & coordinates(xy)[,1] > -4)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Sahara-Mali_")
  }
  
  
  #russian tundra 72 117
aw<-which ( coordinates(xy)[,2] < 73 & coordinates(xy)[,2] > 71 & coordinates(xy)[,1] < 118 & coordinates(xy)[,1] > 116)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Tundra-Rus_")
  }
  
  
 
 
 #kogelberg
 aw<-which ( coordinates(xy)[,2] < -34 & coordinates(xy)[,2] > -35 & coordinates(xy)[,1] < 19 & coordinates(xy)[,1] > 18)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Kogelberg_")
  }
  
#Oslo_fjord
 aw<-which ( coordinates(xy)[,2] < 60 & coordinates(xy)[,2] > 59 & coordinates(xy)[,1] < 10 & coordinates(xy)[,1] > 9)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="Oslo-Fjord-Norway_")
  }
  
#Southern Sweden
 aw<-which ( coordinates(xy)[,2] < 58 & coordinates(xy)[,2] > 57 & coordinates(xy)[,1] < 16 & coordinates(xy)[,1] > 15)
for(pp in 1:10) {
    p<-aw[pp]
    print(p)
    sxy<-round(coordinates(xy)[p,],2)
    stats[p,]<-dopixel(p,PLOT=T,REGION="South-Sweden_")
  }
  
   
  

#  cbind( stats[725895,] ,  stats[1648053,] ,stats[2074751,])
#   
#   }
#   
#  dopixel(1648053,PLOT=T)
#   
#--------------------------------------------------------------------------------------------------
# # create small test dataset to run locally
# xy<-attach.big.matrix("xy.desc")
# cells<-mwhich(xy, cols=1:2, vals=list(c(90,91),c(67,68)), comps=list(c('ge','le'),c('ge','le')), op='AND')
# data<-data[cells,]
# save(data, file="test2.RData")
