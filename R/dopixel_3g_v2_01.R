# New version pulling data from a bigmatrix and writing to a bigmatrix
# first and last leafon and leafoff days are minima/maxima in the sum of the normalised 1st and 2nd derivatives
# see end of script for rationale
# see D:\PhD\NDVI_comparison\GIMMS_NDVI\3g\GIMMS_3g_description.txt for rationale
# on data.e, data.r, and w


#  this is now quite different. 8 July 2015
#  1 make data<0 = 0, it seems like ice and snow data are already set to zero by GIMMS?
#  2 run approx to replace NA data points (data.r==7)
#  3 get a filtered mean on the result of approx
#  4 clean up ends of filter (its circular) by using data rather than filter
#  5 run smooth.spline on filtered data (rather than on approx).
#  this means that NA data are now replaced by the filter, but 0 values

# with Roberts evergreen improvements form 26 August 2015
# with Steve adding some sol radiation information for calculating a GPP index (1 March 2016)

dopixel<-function(p, PLOT=F,REGION="NA") {


  cor.photofac.all<-NA;cor.moistfac.all<-NA;cor.radfac.all<-NA;cor.evi.all<-NA;cv.gpp.all<-NA;cv.evi.all<-NA;
  cor.vpi.photofac.all<-NA; cor.vpi.moistfac.all <-NA; cor.vpi.radfac.all <-NA; cor.vpi.evi.all <-NA;  cv.vpi.all<-NA; 

 
  
  data.e<-floor(data[p,]/10)/1000
  data.r<-data[p,]-floor(data[p,]/10)*10 +1
  
  data.e[data.e<=0]<-NA 
  data.e[which(data.r==4 | data.r==6)]<-0
  data.e[data.r==7]<-NA
    
  NA.length<-length(which(is.na(data.e))) #how many NA data points
  ICE.length<-length(which(data.r==4 | data.r==6))  #how many possible ICE/SNOW data points
  ZERO.length<-length(which(data.e==0))
  
#  FLAG.THIS<-0 #OK
  
#   if( ZERO.length + NA.length >650 ) {stats.ret<-data.frame(matrix(rep(NA, length(names)),nr=1))
#                                       FLAG.THIS<-1  }
if( ZERO.length + NA.length >650 ) {stats.ret<-data.frame(as.list(rep(NA, length(names)))) }
                                        
                                      
   else {
  
  # set a weight for spline fitting, 7 if NA (already assigned NA earlier)
#  w<-c(1,1,0.66,0,0.33,0,0)[data.r]
#  w<-ifelse(is.na(w), 0, w)

  approx.e<-approx(JDAY.x,data.e,xout=JDAY.x,rule=2)$y #use approx to linearly interpolate
  data.e[is.na(data.e)]<-approx.e[is.na(data.e)] #replace NA with linear interpolation
  data.e2<-filter(data.e,filter=rep(1/3,3),side=2,circular=T) #get a filtered running mean:
  data.e2[1]<-data.e[1] #filter is circular so take first value from data
  data.e2[length(data.e)]<-data.e[length(data.e)] #filter is circular so take last value from data
  ss<-smooth.spline(JDAY.x, data.e2, all.knots=T) #smooth spline on filtered data

    lq<-quantile(ss$y, 0.1)
    uq<-quantile(ss$y, 0.9)
    amplitude <- uq-lq
    
    mean.evi<-mean(ss$y) #the mean evi of the time series
    sd.evi<-sd(ss$y) #the standard deviation of the time series
    sum.evi<-sum(ss$y) #sum evi
    
    
    #smoother spline and deriviatives
    sss<-smooth.spline(JDAY.x, approx.e, spar=0.3)
    d1<-predict(sss, sss$x,deriv=1)$y
    d1<-(d1-mean(d1))/sd(d1)
    d2<-predict(sss, sss$x,deriv=2)$y
    d2<-(d2-mean(d2))/sd(d2)
    #lqd<-quantile(d1, 0.05)
    #uqd<-quantile(d1, 0.95)
    
    
    
    ### AVERAGE PEAK/TROUHG/LEAF-ON/LEAF-OFF DAYS
    day.numbers<-sort(unique(JDAY))
    
    peak.day<-day.numbers[by(sss$y,YEAR,which.max)]
    trough.day<-day.numbers[by(sss$y,YEAR,which.min)]
    peak.day<-round( circ.mean(peak.day[-1]))
    trough.day<-round( circ.mean(trough.day[-1]))
    
    #if trough earlier than 1 Jul (data starts 1 Jul 1981), start time series in next year
    trough.x<-ifelse(trough.day>188, trough.day, trough.day+365)
    trough.x<-as.matrix(seq(trough.x, by=365, length.out=33))
    #trough.x<-trough.x + as.numeric(as.Date("1981-01-01"))# add the start date
#    line of code below removed 26.08.2015
#    trough.win.x<-apply(trough.x, 1, FUN=function(x) (x-90):(x+90))
    ###################################################################
    #CHANGED 26.08.2015
    window<-max(16, round((180*amplitude)/2, 0))
    # window of 90 days either side in which to search for annual trough days

    trough.win.x<-apply(trough.x, 1, FUN=function(x) (x-window):(x+window))
    if(class(trough.win.x)=="integer") trough.win.x<-rbind(trough.win.x)
    ####################################################################
    
    td<-NULL; td.x<-NULL; td.evi<-NULL
    
    for(i in 1:33){
      td.win<-which(JDAY.x %in% trough.win.x[,i])
      if(length(td.win)>0){
        trough<-td.win[which.min(sss$y[td.win])]
        td[i]<-JDAY[trough]
        td.x[i]<-JDAY.x[trough]
        td.evi[i]<-ss$y[trough]
      }
      else{ td[i]<-NA; td.x[i]<-NA; td.evi[i]<-NA}
    }
    
    
#---start----STEVE SolRad---------------------------
 # create a spline of SolRad for the entire time series for the pixel
  rad.array.rc<-c( which(lon==txy[p,1]),which(lat==txy[p,2]))
  tmp.array.rc<-c( which(tmp.lon==tmp.txy[p,1]),which(tmp.lat==tmp.txy[p,2]))
  sm.array.rc<-c( which(sm.lon==sm.txy[p,1]),which(sm.lat== -1*sm.txy[p,2])) ## -1 because of weird array packing
  
  SWLAND.array.rc<-c( which(round(land.lon,3)==merra.txy[p,1]),which(round(land.lat,3)==merra.txy[p,2])) #2018 shortwave land
  RZMC.array.rc<-c( which(round(land.lon,3)==merra.txy[p,1]),which(round(land.lat,3)==merra.txy[p,2]))   #2018 rooting zone moisture content
  
  dfss<-NULL
  predict.gpp.all.yrs<-rep(NA,length=length(JDAY.x))
  dftmp<-NULL
  predict.tmp.all.yrs<-rep(NA,length=length(JDAY.x))
  dfsm<-NULL
  predict.sm.all.yrs<-rep(NA,length=length(JDAY.x))
  dfswland<-NULL
  predict.swland.all.yrs<-rep(NA,length=length(JDAY.x))
  dfrzmc<-NULL
  predict.rzmc.all.yrs<-rep(NA,length=length(JDAY.x))

  
  

#  length(which(is.na(dftmp$sr))) > 50
   FLAG.THIS<-1 #set to not OK

  if(length (rad.array.rc) == 2 &  length (sm.array.rc) == 2 &  length (tmp.array.rc) == 2 &
      length (SWLAND.array.rc) == 2 & length (RZMC.array.rc) == 2)   { #if dims <2 FLAT.THIS remains 1 = not OK

   FLAG.THIS<-0  #reset to OK
   dfss<-data.frame(d=m.JDAY.x,sr=rad.array[ rad.array.rc[1],rad.array.rc[2],]) #dim(rad.array) 576 361 420
   dfss<-dfss[complete.cases(dfss),] 
   dftmp<-data.frame(d=tmp.m.JDAY.x,sr=tmp.array[ tmp.array.rc[1],tmp.array.rc[2],]) #dim(rad.array) 
   dftmp<-dftmp[complete.cases(dftmp),] 
   dfsm<-data.frame(d=sm.m.JDAY.x,sr=soilw.array[ sm.array.rc[1],sm.array.rc[2],]) #dim(rad.array) 
   dfsm<-dfsm[complete.cases(dfsm),] 

   dfswland<-data.frame(d=merra.JDAY.x,sr=SWLAND.array[ SWLAND.array.rc[1],SWLAND.array.rc[2],]) 
   dfswland<-dfswland[complete.cases(dfswland),] 
   dfrzmc<-data.frame(d=merra.JDAY.x,sr=RZMC.array[ RZMC.array.rc[1],RZMC.array.rc[2],]) 
   dfrzmc<-dfrzmc[complete.cases(dfrzmc),] 


   
   

   if(dim(dfss)[1]<100 | dim(dftmp)[1] <100 | dim(dfsm)[1] <100 | dim(dfswland)[1] <100 | dim(dfrzmc)[1] <100 ) 
    {FLAG.THIS<-1} # dims too small, so set to not OK  
   else    
    {
    sol.ss<-smooth.spline(dfss$d, dfss$sr, all.knots=T) #smooth spline on sol rad data
    sol.tmp<-smooth.spline(dftmp$d, dftmp$sr, all.knots=T) #smooth spline on tmp data
    sol.sm<-smooth.spline(dfsm$d, dfsm$sr, all.knots=T) #smooth spline on soil moisture data
    sol.swland<-smooth.spline(dfswland$d, dfswland$sr, all.knots=T) #smooth spline on SWLAND
    sol.rzmc<-smooth.spline(dfrzmc$d, dfrzmc$sr, all.knots=T) #smooth spline on rooting zone moisture content
    
    predict.tmp.all.yrs<-predict(sol.tmp,JDAY.x)$y #predict from spline for entire time
    photo.fac.all.yrs<-  -0.242141 + 0.093696*predict.tmp.all.yrs  -0.001768*(predict.tmp.all.yrs)^2  #normalised to max 1
    photo.fac.all.yrs<-pmax(0,photo.fac.all.yrs) #limit this so that photo.fac can not be less than zero
    
    predict.sm.all.yrs<-predict(sol.sm,JDAY.x)$y #predict from spline for entire time
    sm.fac.all.yrs<- pmax(0,pmin(1,(predict.sm.all.yrs-756*0.10)/(756*0.40-756*0.10))) #assumes FC is 40% of max and WP is 10% of max
    sm.fac.all.yrs<-pmax(0,sm.fac.all.yrs) #limit this so that sm.fac can not be less than zero
    
    predict.rzmc.all.yrs<-predict(sol.rzmc,JDAY.x)$y #predict from spline for entire time
    rzmc.fac.all.yrs<- trap1(predict.rzmc.all.yrs,0.1,0.40) #assumes FC is 40% of max and WP is 10% of max
    
    predict.swland.all.yrs<-pmax(0,predict(sol.swland,JDAY.x)$y) #predict from spline for entire time

    predict.sol.all.yrs<-predict(sol.ss,JDAY.x)$y #predict from spline for entire time
    predict.sol.all.yrs<-pmax(0,predict.sol.all.yrs)
    sum.sol.yr.all<-sum(predict.sol.all.yrs)

    
    predict.evi.all.yrs<-pmax(0,predict(ss,JDAY.x)$y) #predict ndvi from spline for entire time
    predict.gpp.all.yrs<-predict.evi.all.yrs * predict.sol.all.yrs * photo.fac.all.yrs * sm.fac.all.yrs 
    predict.vpi.all.yrs<-predict.evi.all.yrs * predict.swland.all.yrs * photo.fac.all.yrs * rzmc.fac.all.yrs 
    



  # 2018 add some correlation calculations
   
#  cor.photofac.all <- cor(predict.gpp.all.yrs,  photo.fac.all.yrs)
#  cor.moistfac.all <- cor(predict.gpp.all.yrs,  sm.fac.all.yrs)
#  cor.radfac.all <- cor(predict.gpp.all.yrs, predict.sol.all.yrs)
#  cor.evi.all <- cor(predict.gpp.all.yrs,  predict.evi.all.yrs)
#  #calling the 2018 index VPI, these are using MERRA Land Solar radiation and Soil Water
#  cor.vpi.photofac.all <- cor(predict.vpi.all.yrs,  photo.fac.all.yrs)
#  cor.vpi.moistfac.all <- cor(predict.vpi.all.yrs,  rzmc.fac.all.yrs)
#  cor.vpi.radfac.all <- cor(predict.vpi.all.yrs, predict.swland.all.yrs)
#  cor.vpi.evi.all <- cor(predict.vpi.all.yrs,  predict.evi.all.yrs)
# 
#  cv.gpp.all <- mycv(predict.gpp.all.yrs)
#  cv.evi.all <- mycv(predict.evi.all.yrs)
#  cv.vpi.all <- mycv(predict.vpi.all.yrs)
 
#  pcr.df<-data.frame(vpi=predict.vpi.all.yrs, temper = photo.fac.all.yrs, moist=rzmc.fac.all.yrs, 
#                      rad =  predict.swland.all.yrs, ndvi =  predict.evi.all.yrs)
# 
#  pcr.vpi<-  pcr(ndvi ~ temper+moist+rad+vpi, ncomp = 4, data = pcr.df,scale=TRUE)
#  
#  cor.vpi.photofac.all <- coef(pcr.vpi)[1]
#  cor.vpi.moistfac.all <- coef(pcr.vpi)[2]
#  cor.vpi.radfac.all <- coef(pcr.vpi)[3]
#  cor.vpi.evi.all <- coef(pcr.vpi)[4]
 
	  
	  
   } # enough data
   } # array indices exist
   
   
   
 
# summary(predict.sol.all.yrs)
# summary(predict.evi.all.yrs)
# summary(photo.fac.all.yrs)
# summary(sm.fac.all.yrs)
# summary(predict.gpp.all.yrs)

  

#---end------STEVE SolRad---------------------------

    
    pd<-NULL; pd.x<-NULL; pd.evi<-NULL;
    elon.m<-NULL; eoff.m<-NULL; elon.m.x<-NULL; eoff.m.x<-NULL; elon.m.evi<-NULL; eoff.m.evi<-NULL
    elon.f<-NULL; eoff.f<-NULL; elon.f.x<-NULL; eoff.f.x<-NULL; elon.f.evi<-NULL; eoff.f.evi<-NULL
    elon.l<-NULL; eoff.l<-NULL; elon.l.x<-NULL; eoff.l.x<-NULL; elon.l.evi<-NULL; eoff.l.evi<-NULL
    sum.evi.yr<-NULL;   
    sum.gpp.yr<-rep(NA,32); #td.gpp<-rep(NA,32); pd.gpp<-rep(NA,32); 
    sum.vpi.yr<-rep(NA,32);
    sum.sol.yr<-rep(NA,32); #STEVE
    td.gpp<-rep(NA,32); pd.gpp<-rep(NA,32); min.gpp<-rep(NA,32); max.gpp<-rep(NA,32);
    amp.gpp<-rep(NA,32); tdg<-rep(NA,32); pdg<-rep(NA,32); tdg.x<-rep(NA,32); pdg.x<-rep(NA,32); #STEVE added sum.gpp.yr
    cor.photofac<-rep(NA,32); cor.moistfac<-rep(NA,32); cor.evi <-rep(NA,32); cor.radfac <-rep(NA,32); cv.gpp <-rep(NA,32);
    cv.evi <-rep(NA,32); 
    cor.vpi.photofac<-rep(NA,32);  cor.vpi.moistfac<-rep(NA,32);  cor.vpi.radfac <-rep(NA,32);  cor.vpi.evi<-rep(NA,32); 
    cv.vpi<- rep(NA,32);    
    
    
    for(i in 1:32){
      
      if(is.na(td[i+1])) {pd[i]<-NA; pd.x[i]<-NA; pd.evi[i]<-NA; elon.m[i]<-NA; elon.m.x[i]<-NA; eoff.m[i]<-NA;
                          eoff.m.x[i]<-NA; elon.f[i]<-NA; elon.f.x[i]<-NA; elon.l[i]<-NA; elon.l.x[i]<-NA;
                          eoff.f[i]<-NA; eoff.f.x[i]<-NA; eoff.l[i]<-NA; eoff.l.x[i]<-NA; sum.evi.yr[i]<-NA;
                          elon.m.evi[i]<-NA; eoff.m.evi[i]<-NA; elon.f.evi[i]<-NA; eoff.f.evi[i]<-NA;
                          elon.l.evi[i]<-NA; eoff.l.evi[i]<-NA;  
                          sum.gpp.yr[i]<-NA; #td.gpp[i]<-NA; pd.gpp[i]<-NA; 
                          sum.vpi.yr[i]<-NA;
                          sum.sol.yr[i]=NA;
                          td.gpp[i]<-NA; pd.gpp[i]<-NA; min.gpp[i]<-NA; max.gpp[i]<-NA;
                          amp.gpp[i]<-NA; tdg[i]<-NA; pdg[i]<-NA; tdg.x[i]<-NA; pdg.x[i]<-NA; #STEVE added sum.gpp.yr
                          cor.photofac[i]<-NA; cor.moistfac[i]<-NA; cor.evi[i]<-NA; cor.radfac[i]<-NA; cv.gpp[i]<-NA; cv.evi[i]<-NA; 
                          cor.vpi.photofac[i]<-NA;  cor.vpi.moistfac[i]<-NA;  cor.vpi.radfac[i] <-NA;  cor.vpi.evi[i]<-NA; 
                          cv.vpi[i] <- NA;
                          } 
      
      else{
        pheno.yr.win<-which(JDAY.x %in% (td.x[i]:(td.x[i+1])-1))
        ###################################################################
        #CHANGED 26.08.2015
        #on-day must fall within the first half of the phenological year
        on.day.win<-pheno.yr.win[1:floor(length(pheno.yr.win)/2)]
        #off-day must fall within the second half of the phenological year
        off.day.win<-(max(pheno.yr.win)-floor(length(pheno.yr.win)/2)):max(pheno.yr.win)
        ###################################################################
#	removed 26.08.2015:        
#         # new, restricting on-day to within 6 months from trough
#         on.day.win<-pheno.yr.win[1:12]
#         # new: off-day within 6 months from next trough
#         off.day.win<-(max(pheno.yr.win)-12):max(pheno.yr.win)
        
        peak<-pheno.yr.win[which.max(sss$y[pheno.yr.win])]
        pd[i]<-JDAY[peak]
        pd.x[i]<-JDAY.x[peak]
        pd.evi[i]<-ss$y[peak]
        
        on.day<-on.day.win[which.max(d1[on.day.win])]
        elon.m[i]<-JDAY[on.day]
        elon.m.x[i]<-JDAY.x[on.day]
        elon.m.evi[i]<-ss$y[on.day]
        
        off.day<-off.day.win[which.min(d1[off.day.win])]
        eoff.m[i]<-JDAY[off.day]
        eoff.m.x[i]<-JDAY.x[off.day]
        eoff.m.evi[i]<-ss$y[off.day]
        
        
        #on.day.f<-pheno.yr.win[min(which(d1[pheno.yr.win]>uqd))]
        on.day.f<-on.day.win[which.max(d1[on.day.win]+d2[on.day.win])]
        elon.f[i]<-JDAY[on.day.f]
        elon.f.x[i]<-JDAY.x[on.day.f]
        elon.f.evi[i]<-ss$y[on.day.f]
        
        #on.day.l<-on.day.win[max(which(d1[on.day.win]>uqd))]
        on.day.l<-on.day.win[which.max(d1[on.day.win]-d2[on.day.win])]
        elon.l[i]<-JDAY[on.day.l]
        elon.l.x[i]<-JDAY.x[on.day.l]
        elon.l.evi[i]<-ss$y[on.day.l]
        
        
        #off.day.f<-pheno.yr.win[min(which(d1[pheno.yr.win]<lqd))]
        off.day.f<-off.day.win[which.min(d1[off.day.win]+d2[off.day.win])]
        eoff.f[i]<-JDAY[off.day.f]
        eoff.f.x[i]<-JDAY.x[off.day.f]
        eoff.f.evi[i]<-ss$y[off.day.f]
        
        #off.day.l<-off.day.win[max(which(d1[off.day.win]<lqd))]
        off.day.l<-off.day.win[which.min(d1[off.day.win]-d2[off.day.win])]
        eoff.l[i]<-JDAY[off.day.l]
        eoff.l.x[i]<-JDAY.x[off.day.l]
        eoff.l.evi[i]<-ss$y[off.day.l]
        
        sum.evi.yr[i]<-sum(ss$y[pheno.yr.win])
        
#---start----STEVE SolRad---------------------------
  if(FLAG.THIS==1) {sum.gpp.yr[i]<-NA;sum.sol.yr[i]<-NA}   else {
  predict.evi.yr<-predict(ss    ,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window
  predict.sol.yr<-predict(sol.ss,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window
  predict.sm.yr<-predict(sol.sm,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window
  predict.tmp.yr<-predict(sol.tmp,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window
  predict.rzmc.yr<-predict(sol.rzmc,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window ADDED-2018
  predict.swland.yr<-pmax(0,predict(sol.swland,(td.x[i]:(td.x[i+1])-1))$y) #predict from spline for time window ADDED-2018
  
  predict.sol.yr<-pmax(0,predict.sol.yr)
  predict.tmp.yr<- -0.242141 +  0.093696*predict.tmp.yr  -0.001768*(predict.tmp.yr)^2  #normalised to max 1
  predict.tmp.yr<- pmax(0,predict.tmp.yr)
  predict.sm.yr<- pmax(0,pmin(1,(predict.sm.yr-756*0.10)/(756*0.40-756*0.10))) #assumes FC is 40% of max and WP is 10% of max
  predict.rzmc.yr<- trap1(predict.rzmc.yr,0.1,0.40) #assumes FC is 40% of max and WP is 10% of max ADDED-2018
  
  
  predict.gpp.yr <- pmax(0,predict.evi.yr) * predict.sol.yr * predict.tmp.yr * predict.sm.yr
  predict.vpi.yr <- pmax(0,predict.evi.yr) * predict.swland.yr * predict.tmp.yr * predict.rzmc.yr #ADDED-2018
  
  sum.gpp.yr[i] <- sum( pmax(0,predict.evi.yr) * predict.sol.yr * predict.tmp.yr * predict.sm.yr) #sum of NDVI * SolarRadiation * TemperatureFac * MoistureFac
  sum.vpi.yr[i] <- sum( pmax(0,predict.evi.yr) * predict.swland.yr * predict.tmp.yr * predict.rzmc.yr) #as above but for 2018 VPI
  sum.sol.yr[i] <- sum( predict.sol.yr)
  
#   #calc VPI using MERRA Land Solrad and Soil Moisture:
#     predict.rzmc.yr<-predict(sol.rzmc,(td.x[i]:(td.x[i+1])-1))$y #predict from spline for time window
#     rzmc.fac.yr<- trap1(predict.rzmc.yr,0.1,0.40) #assumes FC is 40% of max and WP is 10% of max
#     predict.swland.yr<-pmax(0,predict(sol.swland,(td.x[i]:(td.x[i+1])-1))$y) #predict from spline for time window
#     predict.vpi.yr <- pmax(0,predict.evi.yr) * predict.swland.yr * predict.tmp.yr * rzmc.fac.yr
  
#   cor.photofac[i] <- cor(predict.gpp.yr,  predict.tmp.yr)
#   cor.moistfac[i] <- cor(predict.gpp.yr,  predict.sm.yr)
#   cor.radfac[i] <- cor(predict.gpp.yr,  predict.sol.yr)
#   cor.evi[i] <- cor(predict.gpp.yr,  predict.evi.yr)
#   
#   cor.vpi.photofac[i] <- cor(predict.vpi.yr,  predict.tmp.yr)
#   cor.vpi.moistfac[i] <- cor(predict.vpi.yr,  predict.rzmc.yr)
#   cor.vpi.radfac[i] <- cor(predict.vpi.yr,  predict.swland.yr)
#   cor.vpi.evi[i] <- cor(predict.vpi.yr,  predict.evi.yr)
#    
#   cv.gpp[i] <- mycv(predict.gpp.yr)
#   cv.evi[i] <- mycv(predict.evi.yr)
#   cv.vpi[i] <- mycv(predict.vpi.yr)
 
 
# summary(predict.sol.yr)
# summary(predict.evi.yr)
# summary(predict.tmp.yr)
# summary(predict.sm.yr)
# sum.gpp.yr[i]


#this way takes the gpp at peak day  
  pd.gpp[i]<-predict.gpp.all.yrs[peak] #gpp at peak day
  td.gpp[i]<-predict.gpp.all.yrs[trough] #gpp at through day
  
#this way takes the max gpp in pheno.yr.win  
  max.gpp[i]<-max(predict.gpp.all.yrs[pheno.yr.win]) #gpp at peak day
  min.gpp[i]<-min(predict.gpp.all.yrs[pheno.yr.win]) #gpp at through day
  amp.gpp[i]<- max.gpp[i] - min.gpp[i] #amplitude of gpp (using min and max gpp)
  
  peak.gpp   <-pheno.yr.win[which.max(predict.gpp.all.yrs[pheno.yr.win])] #the peak day of gpp
  trough.gpp <-pheno.yr.win[which.min(predict.gpp.all.yrs[pheno.yr.win])] #the trough day  of gpp
  pdg[i]<-JDAY[peak.gpp]   #Julian day of the peak of gpp
  tdg[i]<-JDAY[trough.gpp] #Julian day of the trough of gpp
  pdg.x[i]<-JDAY.x[peak.gpp]   #Julian day of the peak of gpp
  tdg.x[i]<-JDAY.x[trough.gpp] #Julian day of the trough of gpp
  } #end FLAG.THIS
  

 
  
  
#---end------STEVE SolRad---------------------------
        }
    } # end of 32 years
    
#          print(length(sum.evi.yr))
#          print(length(sum.gpp.yr))
    
    amp<-pd.evi-td.evi[-1]
#    amp.gpp<- max.gpp - min.gpp[-1] #STEVE amplitude of gpp (using min and max gpp)
    gsl<-eoff.m.x-elon.m.x
    gsl.long<-eoff.l.x-elon.f.x
    gsl.peak<-eoff.f.x-elon.l.x
    
    
        if(PLOT==TRUE) {
          png(file=paste(REGION,p, "_v8.png", sep=""), width=90,height=10, units="cm", res=300)
          par(mar=c(2,2,4,3))
    
#           col<-ifelse(w==0,1,ifelse(w==0.33,1,ifelse(w==0.66,3,4)))
           col<-1

           plot(JDAY.x,data.e,type="p" ,cex=0.6, pch=16,
               col=col, axes=F, xlab="",ylab="NDVI",
               ylim=c(0,1),
               main=paste("lon", sxy[1], "   lat", sxy[2],
                       "\n c.temp= ", round(cor.photofac.all,2), " c.mois= ", round(cor.moistfac.all,2),
                         " c.rad= ", round(cor.radfac.all,2), " c.ndvi= ", round(cor.evi.all,2), 
                           "\n c.temp= ", round(cor.vpi.photofac.all,2), " c.mois= ", round(cor.vpi.moistfac.all,2),
                         " c.rad= ", round(cor.vpi.radfac.all,2), " c.ndvi= ", round(cor.vpi.evi.all,2) )   
                         
                        )
               

    
          axis(1,las=0,at=seq(0,11671,365),labels=seq(1981,2012,1),cex.axis=1 )
 #         axis(2,at=seq(0,max(data.e), 0.1),labels=seq(0,max(data.e),0.1))
          axis(2,at=seq(0,1, 0.2),labels=seq(0,1,0.2))
          lines(ss$x,ss$y,col="green")
#          lines(sss$x,sss$y)
    
 #         abline(v=seq(0,9481,365), lty=2, col="orange")
          abline(v=td.x, lty=2, col="black")
    
          points(td.x, td.evi, pch=16, col="red", cex=0.6)
          points(elon.f.x,elon.f.evi, col="green", pch=16, cex=0.8)
          points(eoff.l.x,eoff.l.evi, col="red", pch=16,cex=0.8)
          points(elon.l.x,elon.l.evi, col="green", pch=16,cex=0.4)
          points(eoff.f.x,eoff.f.evi, col="red", pch=16,cex=0.4)
          points(elon.m.x,elon.m.evi, col="green", pch=16, cex=0.6)
          points(eoff.m.x,eoff.m.evi, col="red", pch=16, cex=0.6)

         if( FLAG.THIS==0  ) {
         lines(sol.ss$x,sol.ss$y/500,col="orange") #sol rad
         points(dfss$d, dfss$sr/500,pch=16,col="orange",cex=0.6)
         lines(sol.swland$x,sol.swland$y/500,col="orange",lty=2) #new sol rad
         
         lines(JDAY.x,sm.fac.all.yrs,col="blue") #soil moisture factor
         lines(JDAY.x,photo.fac.all.yrs,col="red") # photosynthesis temperature factor
         lines(JDAY.x,predict.gpp.all.yrs/250,col="steelblue") #VPI index all years
         
         lines(JDAY.x,rzmc.fac.all.yrs,col="blue",lty=2) #new soil moisture factor
         lines(JDAY.x,predict.gpp.all.yrs/250,col="steelblue") #VPI index all years
         lines(JDAY.x,predict.vpi.all.yrs/250,col="steelblue",lty=2) #new VPI index all years
         
         
         
         }

        
         
         
         axis(4,at=seq(0,1, length=6),labels=seq(0,500,length=6))

 	  

	  
          legend("topleft", bty="n", ncol=4, inset=c(0,-0.3), xpd=T,
                 c(paste("weight =",seq(0,1,0.33)), paste("leaf-on",1:3),
                   paste("leaf-off",1:3),"trough day", 
                   "ndvi spline", "moisture transp","temperature photo","sol rad","gpp index"),
                 col=c("black","red","green","blue",rep("green",3),rep("red",3),"black",
                      "green","blue","red","orange","steelblue"),
                 lty=c(rep(0,10),2,rep(1,5)), pch=c(rep(1,4),rep(16,6),rep(NA,3)),
                 pt.cex=c(rep(0.4,4),0.8,0.6,0.4,0.4,0.6,0.8,rep(NA,6)))
    
          dev.off()
        } #end plot
    
    devi.ss<-mean(abs(approx.e-ss$y)) # deviation between less smooth spline and data
    devi.sss<-mean(abs(approx.e-sss$y))  # deviation between smoother spline and data
    
    L.JDAY        <-length(JDAY) #length of the input time series
    
    stats.ret<-data.frame(
      lq, uq,  mean.evi, sd.evi, sum.evi,  amplitude,
      peak.day, trough.day,
      devi.ss, devi.sss, L.JDAY,
      NA.length, ICE.length,
      cor.photofac.all, cor.moistfac.all, cor.radfac.all, cor.evi.all, cv.gpp.all, cv.evi.all,
      cor.vpi.photofac.all,cor.vpi.moistfac.all,cor.vpi.radfac.all, cor.vpi.evi.all,cv.vpi.all,  
      as.list(td), as.list(td.x), as.list(td.evi),
      as.list(pd), as.list(pd.x), as.list(pd.evi),
      as.list(elon.m), as.list(elon.m.x), as.list(eoff.m), as.list(eoff.m.x),
      as.list(elon.f), as.list(elon.f.x), as.list(eoff.f), as.list(eoff.f.x),
      as.list(elon.l), as.list(elon.l.x), as.list(eoff.l), as.list(eoff.l.x),
      as.list(sum.evi.yr), as.list(amp), as.list(gsl), as.list(gsl.peak), as.list(gsl.long),
      as.list(elon.f.evi), as.list(elon.m.evi), as.list(elon.l.evi),
      as.list(eoff.f.evi), as.list(eoff.m.evi), as.list(eoff.l.evi),
      as.list(sum.gpp.yr), #as.list(td.gpp),as.list(pd.gpp), 
      as.list(sum.vpi.yr), #as.list(td.gpp),as.list(pd.gpp), 
      as.list(sum.sol.yr),
      as.list(td.gpp), as.list(pd.gpp), as.list(min.gpp), as.list(max.gpp),
      as.list(amp.gpp), 
      as.list(tdg),as.list(pdg),as.list(tdg.x),as.list(pdg.x),
      as.list(cor.photofac), as.list(cor.moistfac),  as.list(cor.radfac), as.list(cor.evi), as.list(cv.gpp), as.list(cv.evi), 
      as.list(cor.vpi.photofac), as.list(cor.vpi.moistfac),  as.list(cor.vpi.radfac), as.list(cor.vpi.evi), as.list(cv.vpi)    
      ) #STEVE added sum.gpp.yr
      
#      print(cbind(sum.gpp.yr,sum.evi.yr,sum.sol.yr))
    
  } #end else ZERO.length
  
  # only needed for testing i.e. when not writing to bigmatrix
  names(stats.ret)<-names
  
  return(stats.ret)
#  print(c(cor.temperature,cor.moistfac,cor.evi))
}



# The first derivative (df1) is max where ndvi increases quickest. The second derivative (df2) is max where df1
# increases the quickest, which should be around where ndvi is at its lowest (i.e. at td). So df2 could be
# useful to indicate the start of the season i.e. where ndvi starts increasing. The difficulty with searching
# for max(df2) within the yearly window is that there might be multiple such maxima within the messy time series.
# Therefore I standardise both df1 and df2, and find the max of their sum. This should fall between where df1
# is max and where df2 is max. Because df1 (in this timeseries) is robust (it correctly finds where ndvi is
# increasin rapdily, because it is calculated from a smoothed timeseries), this method is essentially
# constraining the search for max(df2) to a window around max(df1), because I am maximising their sum.
# Theoretically one should be able to use df2 if the timeseries were smooth enough, i.e. there was only a single
# minimum and maximum in the yearly window. My method has the advantages of being 1) able to use a less-smoothed
# time series (tracking the data more closely) and 2) finding the start of the season somewhere in between
# the point where ndvi starts increasing (max(df2)) and where it inreases most rapidly (max(df1)=elon.m).

# Example code:
# x<-seq(5,11,0.01) #x<-seq(-5,5,0.01)
# y<-sin(x) #y<-1/(1+exp(-x))
# plot(x,y, type="l", lwd=2)
# df1<-predict(smooth.spline(x,y), deriv=1)
# df2<-predict(smooth.spline(x,y), deriv=2)
# lines(df1)
# lines(df2, lty=2)
# df.sum.pos<-df1$y+df2$y
# abline(v=x[which.max(df.sum.pos)], col="green")
# abline(v=x[which.min(df.sum.pos)], lty=2, col="red")
# df.sum.neg<-df1$y-df2$y
# abline(v=x[which.max(df.sum.neg)], lty=2, col="green")
# abline(v=x[which.min(df.sum.neg)], col="red")



#  a<-1:10
#  b<-rnorm(10,0,1)
#  plot(a,b)
#  plot(a,b,main=paste("text",a[1],"text2",round(b[1]), "\ntext3",a[2],"text4",round(b[2])))
