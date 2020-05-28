library(foreign)
library(ncdf4)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(viridis)
library(ggpubr)

wd <- 'D:/Research/globalEmission'

#Common variables, function, and datasets
Qquantiles <- c('Q0', 'Q10', 'Q20', 'Q30', 'Q40', 'Q50', 
                'Q60', 'Q70', 'Q80', 'Q90', 'Q100', 'Qmean', 'Qstd')

hydroBasinID<-read_csv(paste0(wd,'/output/table/hydrobasin/comidHydrobasin4.csv'),
                       col_types='ic')
hydroBasinID<-hydroBasinID[!(hydroBasinID$HYBAS_ID %in% c('0')),]
hydroBasinID_nobasin<-read_csv(paste0(wd,'/output/table/hydrobasin/comidHydrobasin4_nobasin.csv'),
                               col_types='ic')
hydroBasinID_nobasin<-hydroBasinID_nobasin[!(hydroBasinID_nobasin$HYBAS_ID %in% c('0')),]
hydroBasinID<-rbind(hydroBasinID,hydroBasinID_nobasin)
rm(hydroBasinID_nobasin)

hydroBasinSN<-read_csv(paste0(wd,'/output/table/hemis/hydrobasin4.csv'),
                       col_types='cc')
runoff<-read_csv(paste0(wd,'/output/table/hydrobasin/runoffhydrobasin4.csv'),
                 col_types='cn')
runoff<-runoff[!is.na(runoff$runoff),]

hydroBasinFR<-read_csv(paste0(wd,'/output/table/hydrobasin/frclass4.csv'),
                       col_types='cc')
hydroBasinFR_1<-read_csv(paste0(wd,'/output/table/hydrobasin/frclass4_offland.csv'),
                         col_types='cc')
hydroBasinFR<-rbind(hydroBasinFR,hydroBasinFR_1)
rm(hydroBasinFR_1)

#ann q stat extractor
qAnnStatsExt<-function(nc,statsNm){
  qStats<-ncvar_get(nc,statsNm)
  qStats[qStats < 0.000000001]=0.000000001
  qStats<-data.frame(qStats)
  colnames(qStats)<-statsNm
  if (nrow(qStats) == nrow(dbf)){
    return(qStats)
  } else {print('number of rows do not match ...')}
}

#monthly q stats extractor
qStatsExt <- function(nc,statsNm){
  qStats <- ncvar_get(nc, statsNm)
  qStats[qStats < 0.000000001]=0.000000001
  qStats <- data.frame(qStats)
  colnames(qStats) <- paste0(month.abb,statsNm) 
  if (nrow(qStats) == nrow(dbf)){
    return(qStats)
  } else {print('number of rows do not match ...')}
}

for(i in 1:8){
  print(paste0('Processing ',i))
  dbf=read.dbf(paste0(wd,'/dataset/MeritHydroLinPan/pfaf_0',i,'_riv_3sMERIT.dbf'))
  dbf<-dbf[,c('COMID','strmOrder','Slope','Length')]
  
  qAnnStatsNC<-nc_open(paste0(wd,'/dataset/qStat/pfaf_0',i,'_all.nc'))
  dbf<-cbind(dbf,qAnnStatsExt(qAnnStatsNC,Qquantiles[12]))
  names(dbf)[5]<-'yeaQmean'
  
  dbf<-cbind(dbf,qAnnStatsExt(qAnnStatsNC, Qquantiles[13]))
  names(dbf)[6]<-'yeaQstd'
  dbf$yeaQcv<-dbf$yeaQstd/dbf$yeaQmean
  nc_close(qAnnStatsNC)
  
  qStatsNC <- nc_open(paste0(wd, '/dataset/qStat/pfaf_0',i,'.nc'))
  dbf<-cbind(dbf,qStatsExt(qStatsNC, Qquantiles[12]))
  dbf<-cbind(dbf,qStatsExt(qStatsNC, Qquantiles[13]))
  nc_close(qStatsNC)
  
  #convert Qstd to Qrsd
  dbf<-
    dbf%>%mutate(
      JanQcv=JanQstd/JanQmean,
      FebQcv=FebQstd/FebQmean,
      MarQcv=MarQstd/MarQmean,
      AprQcv=AprQstd/AprQmean,
      MayQcv=MayQstd/MayQmean,
      JunQcv=JunQstd/JunQmean,
      JulQcv=JulQstd/JulQmean,
      AugQcv=AugQstd/AugQmean,
      SepQcv=SepQstd/SepQmean,
      OctQcv=OctQstd/OctQmean,
      NovQcv=NovQstd/NovQmean,
      DecQcv=DecQstd/DecQmean
    )
  dbf<-dbf[,!(names(dbf)%in%paste0(month.abb,'Qstd'))]
  uparea<-read_csv(paste0(wd,'/output/table/MeritHydro/uparea_0',i,'.csv'))
  dbf<-left_join(dbf,uparea[,c('COMID','uparea')],by='COMID')
  
  #join HYBAS_ID
  dbf<-left_join(dbf,hydroBasinID,by='COMID')
  if(table(is.na(dbf$HYBAS_ID))[1]!=dim(dbf)[1]) 
  {print(paste(table(is.na(dbf$HYBAS_ID))[2],'Flowlines have no HYBAS_ID'))} else {print('All Flowlines have HYBAS_ID')}
  dbf<-dbf[!is.na(dbf$HYBAS_ID),]
  
  #join S/N
  dbf<-left_join(dbf,hydroBasinSN,by='HYBAS_ID')
  
  #join basin runoff
  dbf<-left_join(dbf,runoff,by='HYBAS_ID')
  print(paste0('There are ',sum(is.na(dbf$runoff)),' flowlines without basin runoff'))
  dbf$Qlevel<-NA
  dbf<-dbf %>% mutate(
    Qlevel=case_when(runoff<50 ~ 'ARID',
                     (runoff<200 & runoff>=50)~'MOD1',
                     (runoff<500 & runoff>=200)~'MOD2',
                     (runoff<1000 & runoff>=500)~'WET1',
                     (runoff>=1000)~'WET2'))
  
  #join frclass
  dbf<-left_join(dbf,hydroBasinFR,by='HYBAS_ID')
  print(paste0('There are ',sum(is.na(dbf$frclass)),' flowlines without FR'))
  dbf[is.na(dbf$frclass),]$COMID
  dbf<-dbf[!is.na(dbf$frclass),]
  dbf<-dbf%>%mutate(wkbasin=paste0(Hemis,Qlevel,substr(frclass,4,6)))
  
  #temp
  temp<-read_csv(paste0(wd,'/output/table/MeritHydro/monTemp_0',i,'.csv'))
  colnames(temp)[2:13]<-paste0(month.abb,'_Ta')
  #calc tw
  temp[,paste0(month.abb,'_Tw')]<-sapply(temp[2:13],function(x){tw=0.67*x+7.45})
  for(col in paste0(month.abb,'_Tw')){temp[temp[,col]<0,col]=0}
  
  #join co2
  co2<-read_csv(paste0(wd,'/output/table/MeritHydro/co2_',str_pad(i,width=2,pad='0'),'_RF_u.csv'))
  co2<-left_join(co2,temp,by='COMID')
  co2[,paste0('co2c_',str_pad(1:12,2,pad='0'))]<-
    mapply(function(x,y){(x-380)*(10^(-1*(-0.00007*y^2 + 0.016*y + 1.11)))},
           co2[,paste0('co2_',str_pad(1:12,2,pad='0'))],co2[,paste0(month.abb,'_Tw')])
  names(co2)[2:13]<-paste0(month.abb,'_pco2')
  names(co2)[38:49]<-paste0(month.abb,'_co2')
  dbf<-left_join(dbf,co2[,c('COMID',paste0(month.abb,'_co2'),paste0(month.abb,'_pco2'))],by='COMID')
  
  #join k
  k<-read_csv(paste0(wd,'/output/table/MeritHydro/k_0',i,'.csv'))
  dbf<-left_join(dbf,k[,c('COMID',paste0(month.abb,'_k'))],by='COMID')
  
  if(i==1){
    df<-dbf
  } else {df<-rbind(df,dbf)}
}
df[df$COMID==61000003,]$uparea<-35
df[df$COMID==31000001,]$uparea<-35
rm(hydroBasinID,runoff,dbf,qStatsNC,i,temp,
   hydroBasinFR,hydroBasinSN,uparea,
   co2,k,col,qAnnStatsExt,qStatsExt,
   qAnnStatsNC,Qquantiles)

df[df$yeaQcv<=quantile(df$yeaQcv,0.001),]$yeaQcv<-quantile(df$yeaQcv,0.001) #0.4,min:0.08
df[df$yeaQcv>=quantile(df$yeaQcv,0.995),]$yeaQcv<-quantile(df$yeaQcv,0.995) #16.7,max:10

#calc annual mean width
df<-
  df%>%mutate(
    yeaWidth=case_when(
      Qlevel%in%c('ARID')~exp(2.1+0.43*log(yeaQmean)),
      Qlevel%in%c('MOD1')~exp(2.1+0.47*log(yeaQmean)),
      Qlevel%in%c('MOD2')~exp(2.24+0.47*log(yeaQmean)),
      Qlevel%in%c('WET1')~exp(2.2+0.45*log(yeaQmean)),
      Qlevel%in%c('WET2')~exp(1.92+0.49*log(yeaQmean))),
    runoffFL=yeaQmean/uparea*3.6*24*365
    )

#calculate WidthExp(exponent of the width-Q relatinship)
#the following procedure prevents unrealistic widthExp estimates
df<-
  df%>%mutate(
    WidthExp=case_when(
      (log(runoffFL)<=1.6&log(yeaQmean)<=(-7)) ~
        0.074*log(yeaQcv)-0.026*1.6-0.011*(-7)+0.36,
      (log(runoffFL)>1.6&log(runoffFL)<=8.5&log(yeaQmean)<= (-7)) ~
        0.074*log(yeaQcv)-0.026*log(runoffFL)-0.011*(-7)+0.36,
      (log(runoffFL)>8.5&log(yeaQmean)<=(-7)) ~
        0.074*log(yeaQcv)-0.026*8.5-0.011*(-7)+0.36,
      (log(runoffFL)<=1.6&log(yeaQmean)>(-7)&log(yeaQmean)<=10) ~
        0.074*log(yeaQcv)-0.026*1.6-0.011*log(yeaQmean)+0.36,
      (log(runoffFL)>1.6&log(runoffFL)<=8.5&log(yeaQmean)>(-7)&log(yeaQmean)<=10) ~
        0.074*log(yeaQcv)-0.026*log(runoffFL)-0.011*log(yeaQmean)+0.36,
      (log(runoffFL)>8.5&log(yeaQmean)>(-7)&log(yeaQmean)<=10) ~
        0.074*log(yeaQcv)-0.026*8.5-0.011*log(yeaQmean)+0.36,
      (log(runoffFL)<=1.6&log(yeaQmean)>10) ~
        0.074*log(yeaQcv)-0.026*1.6-0.011*10+0.36,
      (log(runoffFL)>1.6&log(runoffFL)<=8.5&log(yeaQmean)>10) ~
        0.074*log(yeaQcv)-0.026*log(runoffFL)-0.011*10+0.36,
      (log(runoffFL)>8.5&log(yeaQmean)>10) ~
        0.074*log(yeaQcv)-0.026*8.5-0.011*10+0.36)
  )
df[df$WidthExp<=quantile(df$WidthExp,0.01),]$WidthExp<-quantile(df$WidthExp,0.01)

#calc at-a-station coefs and monthly width
df<-
  df%>%mutate(
    Widthcoef=yeaWidth/(yeaQmean^WidthExp),
    JanWidth=Widthcoef*(JanQmean^WidthExp),
    FebWidth=Widthcoef*(FebQmean^WidthExp),
    MarWidth=Widthcoef*(MarQmean^WidthExp),
    AprWidth=Widthcoef*(AprQmean^WidthExp),
    MayWidth=Widthcoef*(MayQmean^WidthExp),
    JunWidth=Widthcoef*(JunQmean^WidthExp),
    JulWidth=Widthcoef*(JulQmean^WidthExp),
    AugWidth=Widthcoef*(AugQmean^WidthExp),
    SepWidth=Widthcoef*(SepQmean^WidthExp),
    OctWidth=Widthcoef*(OctQmean^WidthExp),
    NovWidth=Widthcoef*(NovQmean^WidthExp),
    DecWidth=Widthcoef*(DecQmean^WidthExp)
  )
df<-df[!df$wkbasin%in%c('NARID15'),] #there is only one flowline, caused by misalignment

#the tibetan quantile
quantile(df[df$HYBAS_ID%in%c('4040838850','4040050470','4040050610','4040050780',
                             '4040050810','4040050910','4040050980','4040051270'),
            paste0(month.abb,'_pco2')]%>%rowMeans(),c(0.05,0.95))

#time dryout
timedryout_calculater<-function(Q){tdryout=1/(1+exp(11+4*log(Q)))}
#calculate timedryout
df[,c('yeatimedryout',paste0(month.abb,'timedryout'))]<-
  sapply(df[,c('yeaQmean',paste0(month.abb,'Qmean'))], timedryout_calculater)
for(mon in month.abb){df[df[,paste0(mon,'Qmean')]>=60,paste0(mon,'timedryout')]=0}

# ####export CO2 efflux####
# df_co2f<-data.frame(df[,c('COMID')])
# names(df_co2f)<-'COMID'
# for(Mon in month.abb){
#   co2Name<-print(paste0(Mon,'_co2'))
#   kName<-print(paste0(Mon,'_k'))
#   co2fName<-print(paste0(Mon,'_co2f'))
#   df_co2f[,co2fName]=df[,co2Name]*df[,kName]*365*12/1000
# }
# #correct extreme values
# for(Mon in month.abb){
#   col<-paste0(Mon,'_co2f')
#   df_co2f[is.na(df_co2f[col]),col]<-median(df_co2f[[col]],na.rm=TRUE)
#   df_co2f[df_co2f[col]>quantile(df_co2f[[col]],0.99,na.rm=TRUE),col]<-quantile(df_co2f[[col]],0.99,na.rm=TRUE)
# }
# write_csv(df_co2f,paste0(wd,'/output/table/MeritHydro/co2f.csv'))

####generate working basins for SO extrapolation####
wkbasins<-count(df,wkbasin)%>%arrange(n)
wkbasins$SOabc<-NA #find out contiguous SO for extrapolation in each basin
for(i in 1:length(wkbasins$wkbasin)){
  print(i)
  df_p<-df[df$wkbasin%in% c(wkbasins[i,]$wkbasin),]
  SOvec<-sort(unique(df_p$strmOrder))
  if(sum(1:5 %in% SOvec)==5){
    wkbasins[i,'SOabc']<-4
  } else if (sum(1:4 %in% SOvec)==4){
    wkbasins[i,'SOabc']<-3
  } else if (sum(1:3 %in% SOvec)==3){
    wkbasins[i,'SOabc']<-2
  } else {wkbasins[i,'SOabc']<-2}
}
rm(df_p)

#join wkbasin prectemp
prectemp<-read_csv(paste0(wd,'/output/table/hydrobasin/basin4TempPrec.csv'))
colnames(prectemp)[2:13]<-paste0(month.abb,'_tavg')
colnames(prectemp)[14:25]<-paste0(month.abb,'_prec')
wkbasins<-left_join(wkbasins,prectemp,by='wkbasin')
rm(prectemp)

####stream order extrapolation and surface area calculation####
for(i in 1:78){
  print(paste0('Processing the ', i, ' basin ...'))
  df_basin<-df[df$wkbasin%in%wkbasins[i,]$wkbasin,]
  for (Mon in month.abb){
    print(Mon)
    widthNM<-paste0(Mon,'Width')
    dryoutNM<-paste0(Mon,'timedryout')
    precNM<-paste0(Mon,'_prec')
    tempNM<-paste0(Mon,'_tavg')
    kNM<-paste0(Mon,'_k')
    co2NM<-paste0(Mon,'_co2')
    pco2NM<-paste0(Mon,'_pco2')
    df_basin_mon<-df_basin[,c('wkbasin','HYBAS_ID','strmOrder','Length',widthNM,dryoutNM,kNM,co2NM,pco2NM)]
    #join prec and temp
    df_basin_mon<-left_join(df_basin_mon,wkbasins[i,c('wkbasin',precNM,tempNM)],by='wkbasin')
    colnames(df_basin_mon)[5]<-'Width'
    colnames(df_basin_mon)[6]<-'timedryout'
    colnames(df_basin_mon)[7]<-'k'
    colnames(df_basin_mon)[8]<-'co2'
    colnames(df_basin_mon)[9]<-'pco2'
    colnames(df_basin_mon)[10]<-'prec'
    colnames(df_basin_mon)[11]<-'temp'
    
    df_basin_mon<-df_basin_mon%>%mutate(
      isChannel=as.numeric(df_basin_mon$Width>=0.3),
      surfArea=Length*Width*isChannel/1000000,
      ephemArea=timedryout*Width*Length*isChannel/1000000
    )
    
    # df_basin_mon[is.na(df_basin_mon$co2E),]%>%View()
    ####use ray2013 method to predict ephemeralness
    #1
    df_basin_mon<-
      df_basin_mon%>%mutate(
        percInterm_ray13=case_when(strmOrder>4~0,
                                   strmOrder==4~((-0.005)*prec+0.023*temp+0.27),
                                   strmOrder==3~((-0.008)*prec+0.028*temp+0.44),
                                   strmOrder==2~((-0.009)*prec+0.029*temp+0.61),
                                   strmOrder==1~((-0.009)*prec+0.029*temp+0.77)),
        timedryout_ray13=case_when(strmOrder>4~0,
                                   strmOrder==4~((-0.0028)*prec+0.012*temp+0.127),
                                   strmOrder==3~((-0.0018)*prec+0.011*temp+0.058),
                                   strmOrder==2~((-0.0021)*prec+0.011*temp+0.088),
                                   strmOrder==1~((-0.0019)*prec+0.017*temp+0.026)))
    #2
    df_basin_mon<-
      df_basin_mon%>%mutate(
        percInterm_ray13=case_when(percInterm_ray13>0.9~0.9,
                                   (percInterm_ray13<=0.9)&(percInterm_ray13>=0)~percInterm_ray13,
                                   percInterm_ray13<0~0),
        timedryout_ray13=case_when(timedryout_ray13>0.9~0.9,
                                   (timedryout_ray13<=0.9)&(timedryout_ray13>=0)~timedryout_ray13,
                                   timedryout_ray13<0~0),
        ephemArea_ray13=percInterm_ray13*timedryout_ray13*Width*Length*isChannel/1000000
      )
  
    #SO summary
    SOsum<-
      df_basin_mon%>%group_by(strmOrder)%>%summarise(
        length_km=sum(Length*isChannel/1000), #km
        width_m=mean(Width*isChannel),
        area_km2=sum(surfArea),
        ephemArea_km2=sum(ephemArea),
        ephemAreaRatio=ephemArea_km2/area_km2,
        ephemArea_km2_ray13=sum(ephemArea_ray13),
        ephemAreaRatio_ray13=ephemArea_km2_ray13/area_km2,
        k=sum(k*(surfArea-ephemArea),na.rm=TRUE)/area_km2,
        co2=sum(co2*(surfArea-ephemArea),na.rm=TRUE)/area_km2,
        pco2=sum(pco2*(surfArea-ephemArea),na.rm=TRUE)/area_km2)
    
    totEphemArea_1<-sum(SOsum$ephemArea_km2)
    totEphemArea_1_ray13<-sum(SOsum$ephemArea_km2_ray13)
    
    #extrapolation
    #use width extrapolation to find out endSO
    #predict width for non-captured SO
    fit<-lm(log(width_m)~strmOrder,data=SOsum[1:wkbasins[i,]$SOabc,])
    fitsum<-summary(fit)
    r2_width<-fitsum$r.squared
    intercept<-round(fit$coefficients[1],5)
    slope<-round(fit$coefficients[2],5)
    endSO<-(log(0.3)-intercept)/slope
    endSOi<-0:ceiling(endSO)
    SO<-c(endSOi,endSO)
    # SO<-c(endSOi,floor(endSO))
    width_extrap<-data.frame(strmOrder=SO)
    width_extrap$width_m<-exp(predict(fit,width_extrap))
    
    #predict totLength for non-captured SOs
    #using endSO estimated from width extrapolation and length scaling relationship across SO
    fit <- lm(log(length_km)~strmOrder,data=SOsum[1:wkbasins[i,]$SOabc,])
    fitsum<-summary(fit)
    r2_length<-fitsum$r.squared
    length_extrap<-data.frame(strmOrder=SO)
    length_extrap$length_km<-exp(predict(fit,length_extrap))
    
    #extrapolate k
    fit<-lm(k~strmOrder,data=SOsum[1:wkbasins[i,]$SOabc,])
    fitsum<-summary(fit)
    r2_k<-fitsum$r.squared
    k_extrap<-data.frame(strmOrder=SO)
    if(fitsum$coefficients[2]>0){k_extrap$k<-mean(SOsum$k,na.rm=TRUE)}else{k_extrap$k<-predict(fit,k_extrap)}
    
    # find the largest non-zero SO for ephemeralness
    SOabc_e<-wkbasins[i,]$SOabc
    while((SOsum[SOsum$strmOrder %in% c(SOabc_e),]$ephemAreaRatio==0)){
      SOabc_e<-SOabc_e-1
      if(SOabc_e<1){break}}
    
    if(SOabc_e>=2){
      fit <- lm(ephemAreaRatio~strmOrder,data=SOsum[1:SOabc_e,])
      fitsum<-summary(fit)
      r2_ephemArea<-fitsum$r.squared
      emphemAreaRatio_extrap<-data.frame(strmOrder=SO)
      emphemAreaRatio_extrap$ephemAreaRatio<-predict(fit,emphemAreaRatio_extrap)
      totEphemArea_2<-
        sum(emphemAreaRatio_extrap$ephemAreaRatio*width_extrap$width_m*length_extrap$length_km/1000)
    } else{
      r2_ephemArea<-NA
      totEphemArea_2<-NA
      emphemAreaRatio_extrap<-data.frame(strmOrder=SO)
      emphemAreaRatio_extrap$ephemAreaRatio<-0}
    
    if(SOabc_e>=2){
      fit <- lm(ephemAreaRatio_ray13~strmOrder,data=SOsum[1:SOabc_e,])
      fitsum<-summary(fit)
      r2_ephemArea_ray13<-fitsum$r.squared
      if(is.na(r2_ephemArea_ray13)){
        r2_ephemArea_ray13='No Ephem'
        emphemAreaRatio_ray13_extrap<-data.frame(strmOrder=SO)
        emphemAreaRatio_ray13_extrap$ephemAreaRatio_ray13<-0
        totEphemArea_2_ray13=0
      } else{
        emphemAreaRatio_ray13_extrap<-data.frame(strmOrder=SO)
        emphemAreaRatio_ray13_extrap$ephemAreaRatio_ray13<-predict(fit,emphemAreaRatio_ray13_extrap)
        totEphemArea_2_ray13<-
          sum(emphemAreaRatio_ray13_extrap$ephemAreaRatio_ray13*width_extrap$width_m*length_extrap$length_km/1000)
      }
    } else{
      r2_ephemArea_ray13<-NA
      totEphemArea_2_ray13<-NA}
    
    extrap<-cbind(wkbasins[i,]$wkbasin,width_extrap,length_extrap[,2],emphemAreaRatio_extrap[,2],emphemAreaRatio_ray13_extrap[,2],k_extrap[,2])
    names(extrap)<-c('wkbasin','strmOrder','width_m','length_km','ephemAreaRatio','ephemAreaRatio_ray13','k')
    extrap$wkbasin<-as.character(extrap$wkbasin)
    extrap['co2']<-SOsum[SOsum$strmOrder==1,'co2']
    extrap['pco2']<-SOsum[SOsum$strmOrder==1,'pco2']
    extrap<-extrap%>%mutate(
      area_km2=width_m*length_km/1000,
      ephemArea_km2=area_km2*ephemAreaRatio,
      ephemArea_km2_ray13=area_km2*ephemAreaRatio_ray13
    )
    #correct for negative effective area
    extrap[extrap$ephemArea_km2>extrap$area_km2,]$ephemArea_km2<-extrap[extrap$ephemArea_km2>extrap$area_km2,]$area_km2
    extrap[extrap$ephemArea_km2_ray13>extrap$area_km2,]$ephemArea_km2_ray13<-extrap[extrap$ephemArea_km2_ray13>extrap$area_km2,]$area_km2
    extrap<-extrap%>%mutate(
      effecArea_km2=area_km2-ephemArea_km2,
      co2F=co2*k*12*0.365,#gC/m2/yr
      co2E=co2*k*effecArea_km2*365*12/1000000#GgC/yr
    )
    
    totLength_2<-sum(extrap$length_km) #total extrapolated flowline length
    totArea_2<-sum(extrap$area_km2)# total extrapolated surface area
    totEphemArea_2<-sum(extrap$ephemArea_km2)
    totEffectArea_2<-sum(extrap$effecArea_km2)
    totEphemArea_2_ray13<-sum(extrap$ephemArea_km2_ray13)
    if(sum(extrap$effecArea_km2!=0)){k_2<-sum(extrap$k*extrap$effecArea_km2)/sum(extrap$effecArea_km2)}else{k_2=0}
    co2_2<-mean(extrap$co2)
    pco2_2<-mean(extrap$pco2)
    co2F_2<-if(sum(extrap$effecArea_km2!=0)){sum(extrap$co2F*extrap$effecArea_km2)/totEffectArea_2}else{co2F_2=0}
    co2E_2<-sum(extrap$co2E)
    
    basinMonArea<-data.frame(basinCode=wkbasins[i,]$wkbasin,Mon,
                             endSO,r2_width,r2_length,r2_ephemArea,r2_ephemArea_ray13,r2_k,totEphemArea_1,
                             totEphemArea_1_ray13,
                             totLength_2,totArea_2,totEphemArea_2,totEphemArea_2_ray13,totEffectArea_2,
                             k_2,co2_2,pco2_2,co2F_2,co2E_2)
    
    #combining results
    if(i==1 & Mon=='Jan'){basinArea<-basinMonArea} else {basinArea<-rbind(basinArea,basinMonArea)}
  }
}
rm(basinMonArea, df_basin,df_basin_mon,emphemAreaRatio_extrap,emphemAreaRatio_ray13_extrap,extrap,fit,fitsum,k_extrap,length_extrap,
   SOsum,width_extrap,co2_2,co2E_2,co2NM,dryoutNM,endSO,endSOi,i,intercept,k_2,kNM,Mon,precNM,totEphemArea_1_ray13,
   r2_ephemArea,r2_ephemArea_ray13,r2_k,r2_length,r2_width,slope,SO,SOabc_e,tempNM,totArea_2,totEphemArea_1,
   totEffectArea_2, totEphemArea_2,totEphemArea_2_ray13,totLength_2,widthNM,co2F_2,pco2NM,pco2_2)
# write_csv(basinArea,paste0(wd,'/output/table/regionSurfArea/basinArea.csv'))



EphemAreaSum<-
  basinArea%>%group_by(Mon)%>%summarise(
    ephemArea=sum(totEphemArea_1)+sum(totEphemArea_2),
    ephemArea_ray13=sum(totEphemArea_1_ray13)+sum(totEphemArea_2_ray13))
write_csv(EphemAreaSum,paste0(wd,'/output/table/regionSurfArea/EphemArea.csv'))


####replacing yeaWidth with GRWL width where available####
GRWLwidth<-read_csv(paste0(wd,'/output/table/GRWL/GRWLwidthHydroBASIN4_30mplus.csv'))
GRWLwidth<-GRWLwidth[,c('COMID','width_mean')]
GRWLwidth<-GRWLwidth%>%group_by(COMID)%>%summarise(width_mean=min(width_mean))
GRWLwidth<-GRWLwidth[GRWLwidth$width_mean>=90,]
nrow(GRWLwidth)/nrow(df) #8.5%,4.3%

df_1<-df[!df$COMID%in%GRWLwidth$COMID,]
df_2<-df[df$COMID%in%GRWLwidth$COMID,]
rm(df)
gc()

df_2<-left_join(df_2,GRWLwidth,by='COMID')
df_2$rt<-df_2$width_mean/df_2$yeaWidth
if(sum(df_2$rt>2)){df_2[df_2$rt>2,]$rt=2}
# rtm<-median(df_2[df_2$rt<3,]$rt,na.rm=TRUE) #ratio(rt) too high (e.g.>3) is because misalignment.
# if(is.na(rtm)){rtm=1}
df_2<-
  df_2%>%mutate(
    JanWidth=JanWidth*rt,
    FebWidth=FebWidth*rt,
    MarWidth=MarWidth*rt,
    AprWidth=AprWidth*rt,
    MayWidth=MayWidth*rt,
    JunWidth=JunWidth*rt,
    JulWidth=JulWidth*rt,
    AugWidth=AugWidth*rt,
    SepWidth=SepWidth*rt,
    OctWidth=OctWidth*rt,
    NovWidth=NovWidth*rt,
    DecWidth=DecWidth*rt)
df_2<-df_2[,!names(df_2)%in%c('yeaWidth')]
names(df_2)[names(df_2)%in%c('width_mean')]<-'yeaWidth'
df_2<-df_2[,names(df_1)]
if(sum(names(df_1)==names(df_2))==ncol(df_1)){
  df<-rbind(df_1,df_2)}else{print('GRWL width is not added sucessfully!')}
rm(df_1,df_2)
gc()

####join hydroBAS atts####
df<-df[!(df$HYBAS_ID%in%c("1040040050","2040059370","5040055420")),]#basins have no valid rivArea
df<-df[!is.na(df$Apr_k),]
#linking wkbasin to HydroSHEDS04 basins
hydroBAS<-df%>%group_by(HYBAS_ID)%>%summarise(wkbasin=wkbasin[1])
names(hydroBAS)[2]<-'basinCode'
# write.csv(hydroBAS,paste0(wd,'/output/table/flowregime/hydroBAS.csv'))
#join basinCentroid
basinCentroid<-read_csv(paste0(wd,'/output/table/hydrobasin/hydrobasin4_centroid.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
hydroBAS<-left_join(hydroBAS,basinCentroid,by='HYBAS_ID')
#give climate zones
hydroBAS<-
  hydroBAS%>%mutate(climzone=case_when((Lat>56)~'Polar',
                                       (Lat<=56&Lat>23.5)~'North Temperate',
                                       (Lat<=23.5&Lat>=-23.5)~'Tropical',
                                       (Lat<=-23.5)~'South Temperate'))
hydroBAS$climzone<-factor(hydroBAS$climzone,levels=c('Polar','North Temperate','Tropical','South Temperate'))
#join basin area
basin04area<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/area04.csv'),col_types='cd')#km2
names(basin04area)[2]<-'basinArea'
hydroBAS<-left_join(hydroBAS,basin04area,by='HYBAS_ID')
rm(basin04area)
#sum up basin areas belonging to the same WKbasin
WKbasinArea<-hydroBAS%>%group_by(basinCode)%>%summarise(wkbasinArea=sum(basinArea))
hydroBAS<-left_join(hydroBAS,WKbasinArea,by="basinCode")
hydroBAS$areaRatio<-hydroBAS$basinArea/hydroBAS$wkbasinArea
rm(WKbasinArea)
#join prectemp
prectemp<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/tempPrep.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
prectemp<-prectemp[,c('HYBAS_ID','annTemp','annPrec')] #degree celcius, mm/yr
hydroBAS<-left_join(hydroBAS,prectemp,by='HYBAS_ID')
rm(prectemp)
#join runoff
runoff<-read_csv(paste0(wd,'/output/table/hydrobasin/runoffhydrobasin4.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
runoff<-runoff%>%mutate(
  wetness=case_when(runoff<50~'Arid',
                    runoff>=50&runoff<500~'Mod',
                    runoff>=500~'Wet')
)
hydroBAS<-left_join(hydroBAS,runoff,by='HYBAS_ID')
rm(runoff)
#join ann gpp npp
anngppnpp<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/annPP.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
hydroBAS<-left_join(hydroBAS,anngppnpp,by='HYBAS_ID')
rm(anngppnpp)
#join mon gpp npp
mongppnpp<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/monPP_new.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
hydroBAS<-left_join(hydroBAS,mongppnpp,by='HYBAS_ID')
rm(mongppnpp)
#join soil respiration
soilResp<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/soilResp_buffer.csv'),col_types=cols(.default='d',HYBAS_ID='c'))
#gCm-2yr-1
soilResp[,paste0('pRS_',str_pad(1:12,2,pad='0'))]<-sapply(soilResp[,paste0('pRS_',str_pad(1:12,2,pad='0'))],function(x){x*365})
hydroBAS<-left_join(hydroBAS,soilResp,by='HYBAS_ID')
rm(soilResp)
####join part1####
#join pc02_1
pco2_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(pco2_Jan=sum(Jan_pco2*Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3))/
              sum(Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3)),
            pco2_Feb=sum(Feb_pco2*Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3))/
              sum(Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3)),
            pco2_Mar=sum(Mar_pco2*Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3))/
              sum(Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3)),
            pco2_Apr=sum(Apr_pco2*Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3))/
              sum(Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3)),
            pco2_May=sum(May_pco2*Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3))/
              sum(Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3)),
            pco2_Jun=sum(Jun_pco2*Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3))/
              sum(Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3)),
            pco2_Jul=sum(Jul_pco2*Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3))/
              sum(Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3)),
            pco2_Aug=sum(Aug_pco2*Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3))/
              sum(Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3)),
            pco2_Sep=sum(Sep_pco2*Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3))/
              sum(Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3)),
            pco2_Oct=sum(Oct_pco2*Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3))/
              sum(Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3)),
            pco2_Nov=sum(Nov_pco2*Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3))/
              sum(Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3)),
            pco2_Dec=sum(Dec_pco2*Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3))/
              sum(Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3)))
hydroBAS_res1<-left_join(hydroBAS,pco2_1,by='HYBAS_ID')
rm(pco2_1)
gc()

#k_1
k_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(k_Jan=sum(Jan_k*Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3))/
              sum(Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3)),
            k_Feb=sum(Feb_k*Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3))/
              sum(Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3)),
            k_Mar=sum(Mar_k*Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3))/
              sum(Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3)),
            k_Apr=sum(Apr_k*Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3))/
              sum(Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3)),
            k_May=sum(May_k*Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3))/
              sum(Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3)),
            k_Jun=sum(Jun_k*Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3))/
              sum(Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3)),
            k_Jul=sum(Jul_k*Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3))/
              sum(Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3)),
            k_Aug=sum(Aug_k*Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3))/
              sum(Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3)),
            k_Sep=sum(Sep_k*Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3))/
              sum(Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3)),
            k_Oct=sum(Oct_k*Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3))/
              sum(Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3)),
            k_Nov=sum(Nov_k*Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3))/
              sum(Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3)),
            k_Dec=sum(Dec_k*Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3))/
              sum(Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3)))
hydroBAS_res1<-left_join(hydroBAS_res1,k_1,by='HYBAS_ID')
rm(k_1)
gc()

#co2F_1
co2F_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(co2F_Jan=sum(Jan_co2*Jan_k*365*12/1000*Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3))/
              sum(Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3)),
            co2F_Feb=sum(Feb_co2*Feb_k*365*12/1000*Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3))/
              sum(Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3)),
            co2F_Mar=sum(Mar_co2*Mar_k*365*12/1000*Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3))/
              sum(Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3)),
            co2F_Apr=sum(Apr_co2*Apr_k*365*12/1000*Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3))/
              sum(Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3)),
            co2F_May=sum(May_co2*May_k*365*12/1000*Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3))/
              sum(Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3)),
            co2F_Jun=sum(Jun_co2*Jun_k*365*12/1000*Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3))/
              sum(Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3)),
            co2F_Jul=sum(Jul_co2*Jul_k*365*12/1000*Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3))/
              sum(Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3)),
            co2F_Aug=sum(Aug_co2*Aug_k*365*12/1000*Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3))/
              sum(Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3)),
            co2F_Sep=sum(Sep_co2*Sep_k*365*12/1000*Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3))/
              sum(Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3)),
            co2F_Oct=sum(Oct_co2*Oct_k*365*12/1000*Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3))/
              sum(Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3)),
            co2F_Nov=sum(Nov_co2*Nov_k*365*12/1000*Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3))/
              sum(Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3)),
            co2F_Dec=sum(Dec_co2*Dec_k*365*12/1000*Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3))/
              sum(Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3)))
hydroBAS_res1<-left_join(hydroBAS_res1,co2F_1,by='HYBAS_ID')
rm(co2F_1)
gc()

#join totArea_1
totArea_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(rivArea_Jan=sum(Length*JanWidth/1000000*(JanWidth>=0.3)),
            rivArea_Feb=sum(Length*FebWidth/1000000*(FebWidth>=0.3)),
            rivArea_Mar=sum(Length*MarWidth/1000000*(MarWidth>=0.3)),
            rivArea_Apr=sum(Length*AprWidth/1000000*(AprWidth>=0.3)),
            rivArea_May=sum(Length*MayWidth/1000000*(MayWidth>=0.3)),
            rivArea_Jun=sum(Length*JunWidth/1000000*(JunWidth>=0.3)),
            rivArea_Jul=sum(Length*JulWidth/1000000*(JulWidth>=0.3)),
            rivArea_Aug=sum(Length*AugWidth/1000000*(AugWidth>=0.3)),
            rivArea_Sep=sum(Length*SepWidth/1000000*(SepWidth>=0.3)),
            rivArea_Oct=sum(Length*OctWidth/1000000*(OctWidth>=0.3)),
            rivArea_Nov=sum(Length*NovWidth/1000000*(NovWidth>=0.3)),
            rivArea_Dec=sum(Length*DecWidth/1000000*(DecWidth>=0.3)))
hydroBAS_res1<-left_join(hydroBAS_res1,totArea_1,by='HYBAS_ID')
rm(totArea_1)
gc()

#Join ephemArea_1
ephemArea_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(ephemArea_Jan=sum(Length*JanWidth/1000000*(Jantimedryout)*(JanWidth>=0.3)),
            ephemArea_Feb=sum(Length*FebWidth/1000000*(Febtimedryout)*(FebWidth>=0.3)),
            ephemArea_Mar=sum(Length*MarWidth/1000000*(Martimedryout)*(MarWidth>=0.3)),
            ephemArea_Apr=sum(Length*AprWidth/1000000*(Aprtimedryout)*(AprWidth>=0.3)),
            ephemArea_May=sum(Length*MayWidth/1000000*(Maytimedryout)*(MayWidth>=0.3)),
            ephemArea_Jun=sum(Length*JunWidth/1000000*(Juntimedryout)*(JunWidth>=0.3)),
            ephemArea_Jul=sum(Length*JulWidth/1000000*(Jultimedryout)*(JulWidth>=0.3)),
            ephemArea_Aug=sum(Length*AugWidth/1000000*(Augtimedryout)*(AugWidth>=0.3)),
            ephemArea_Sep=sum(Length*SepWidth/1000000*(Septimedryout)*(SepWidth>=0.3)),
            ephemArea_Oct=sum(Length*OctWidth/1000000*(Octtimedryout)*(OctWidth>=0.3)),
            ephemArea_Nov=sum(Length*NovWidth/1000000*(Novtimedryout)*(NovWidth>=0.3)),
            ephemArea_Dec=sum(Length*DecWidth/1000000*(Dectimedryout)*(DecWidth>=0.3)))
hydroBAS_res1<-left_join(hydroBAS_res1,ephemArea_1,by='HYBAS_ID')
rm(ephemArea_1)
gc()

#EffectArea_1
effectArea_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(effectArea_Jan=sum(Length*JanWidth/1000000*(1-Jantimedryout)*(JanWidth>=0.3)),
            effectArea_Feb=sum(Length*FebWidth/1000000*(1-Febtimedryout)*(FebWidth>=0.3)),
            effectArea_Mar=sum(Length*MarWidth/1000000*(1-Martimedryout)*(MarWidth>=0.3)),
            effectArea_Apr=sum(Length*AprWidth/1000000*(1-Aprtimedryout)*(AprWidth>=0.3)),
            effectArea_May=sum(Length*MayWidth/1000000*(1-Maytimedryout)*(MayWidth>=0.3)),
            effectArea_Jun=sum(Length*JunWidth/1000000*(1-Juntimedryout)*(JunWidth>=0.3)),
            effectArea_Jul=sum(Length*JulWidth/1000000*(1-Jultimedryout)*(JulWidth>=0.3)),
            effectArea_Aug=sum(Length*AugWidth/1000000*(1-Augtimedryout)*(AugWidth>=0.3)),
            effectArea_Sep=sum(Length*SepWidth/1000000*(1-Septimedryout)*(SepWidth>=0.3)),
            effectArea_Oct=sum(Length*OctWidth/1000000*(1-Octtimedryout)*(OctWidth>=0.3)),
            effectArea_Nov=sum(Length*NovWidth/1000000*(1-Novtimedryout)*(NovWidth>=0.3)),
            effectArea_Dec=sum(Length*DecWidth/1000000*(1-Dectimedryout)*(DecWidth>=0.3)))
#get ice coverage
icecov<-read_csv(paste0(wd,'/output/table/HydroBASINSatts/iceout.csv'),
                 col_types=cols(.default='d',HYBAS_ID='c'))
names(icecov)<-c('HYBAS_ID',paste0('iceCov_',month.abb))
effectArea_1<-left_join(effectArea_1,icecov,by='HYBAS_ID')
#calculate iceCovered area
effectArea_1<-
  effectArea_1%>%mutate(
  icecovArea_Jan=effectArea_Jan*iceCov_Jan,
  icecovArea_Feb=effectArea_Feb*iceCov_Feb,
  icecovArea_Mar=effectArea_Mar*iceCov_Mar,
  icecovArea_Apr=effectArea_Apr*iceCov_Apr,
  icecovArea_May=effectArea_May*iceCov_May,
  icecovArea_Jun=effectArea_Jun*iceCov_Jun,
  icecovArea_Jul=effectArea_Jul*iceCov_Jul,
  icecovArea_Aug=effectArea_Aug*iceCov_Aug,
  icecovArea_Sep=effectArea_Sep*iceCov_Sep,
  icecovArea_Oct=effectArea_Oct*iceCov_Oct,
  icecovArea_Nov=effectArea_Nov*iceCov_Nov,
  icecovArea_Dec=effectArea_Dec*iceCov_Dec,
  effectArea_Jan=effectArea_Jan-icecovArea_Jan,
  effectArea_Feb=effectArea_Feb-icecovArea_Feb,
  effectArea_Mar=effectArea_Mar-icecovArea_Mar,
  effectArea_Apr=effectArea_Apr-icecovArea_Apr,
  effectArea_May=effectArea_May-icecovArea_May,
  effectArea_Jun=effectArea_Jun-icecovArea_Jun,
  effectArea_Jul=effectArea_Jul-icecovArea_Jul,
  effectArea_Aug=effectArea_Aug-icecovArea_Aug,
  effectArea_Sep=effectArea_Sep-icecovArea_Sep,
  effectArea_Oct=effectArea_Oct-icecovArea_Oct,
  effectArea_Nov=effectArea_Nov-icecovArea_Nov,
  effectArea_Dec=effectArea_Dec-icecovArea_Dec
)
#join iceCovArea and effectArea
hydroBAS_res1<-left_join(hydroBAS_res1,effectArea_1[,c('HYBAS_ID',paste0('icecovArea_',month.abb),
                                    paste0('effectArea_',month.abb))],by='HYBAS_ID')

rm(effectArea_1)
gc()

#join co2E_1,unit: 10^9gCyr-1
co2E_1<-
  df%>%group_by(HYBAS_ID)%>%
  summarise(co2E_Jan=sum(Jan_co2*Jan_k*365*12/1000*Length*JanWidth/1000000000*(1-Jantimedryout)*(JanWidth>=0.3)),
            co2E_Feb=sum(Feb_co2*Feb_k*365*12/1000*Length*FebWidth/1000000000*(1-Febtimedryout)*(FebWidth>=0.3)),
            co2E_Mar=sum(Mar_co2*Mar_k*365*12/1000*Length*MarWidth/1000000000*(1-Martimedryout)*(MarWidth>=0.3)),
            co2E_Apr=sum(Apr_co2*Apr_k*365*12/1000*Length*AprWidth/1000000000*(1-Aprtimedryout)*(AprWidth>=0.3)),
            co2E_May=sum(May_co2*May_k*365*12/1000*Length*MayWidth/1000000000*(1-Maytimedryout)*(MayWidth>=0.3)),
            co2E_Jun=sum(Jun_co2*Jun_k*365*12/1000*Length*JunWidth/1000000000*(1-Juntimedryout)*(JunWidth>=0.3)),
            co2E_Jul=sum(Jul_co2*Jul_k*365*12/1000*Length*JulWidth/1000000000*(1-Jultimedryout)*(JulWidth>=0.3)),
            co2E_Aug=sum(Aug_co2*Aug_k*365*12/1000*Length*AugWidth/1000000000*(1-Augtimedryout)*(AugWidth>=0.3)),
            co2E_Sep=sum(Sep_co2*Sep_k*365*12/1000*Length*SepWidth/1000000000*(1-Septimedryout)*(SepWidth>=0.3)),
            co2E_Oct=sum(Oct_co2*Oct_k*365*12/1000*Length*OctWidth/1000000000*(1-Octtimedryout)*(OctWidth>=0.3)),
            co2E_Nov=sum(Nov_co2*Nov_k*365*12/1000*Length*NovWidth/1000000000*(1-Novtimedryout)*(NovWidth>=0.3)),
            co2E_Dec=sum(Dec_co2*Dec_k*365*12/1000*Length*DecWidth/1000000000*(1-Dectimedryout)*(DecWidth>=0.3)))
#correct for ice cover
co2E_1<-left_join(co2E_1,icecov,by='HYBAS_ID')
co2E_1<-co2E_1%>%mutate(
  co2E_Jan=co2E_Jan*(1-iceCov_Jan),
  co2E_Feb=co2E_Feb*(1-iceCov_Feb),
  co2E_Mar=co2E_Mar*(1-iceCov_Mar),
  co2E_Apr=co2E_Apr*(1-iceCov_Apr),
  co2E_May=co2E_May*(1-iceCov_May),
  co2E_Jun=co2E_Jun*(1-iceCov_Jun),
  co2E_Jul=co2E_Jul*(1-iceCov_Jul),
  co2E_Aug=co2E_Aug*(1-iceCov_Aug),
  co2E_Sep=co2E_Sep*(1-iceCov_Sep),
  co2E_Oct=co2E_Oct*(1-iceCov_Oct),
  co2E_Nov=co2E_Nov*(1-iceCov_Nov),
  co2E_Dec=co2E_Dec*(1-iceCov_Dec)
)
hydroBAS_res1<-left_join(hydroBAS_res1,co2E_1[,1:13],by='HYBAS_ID')
rm(co2E_1)
gc()

####join part2####
#join res2
cols=c('pco2_2','k_2','co2F_2','totArea_2','totEphemArea_2','totEffectArea_2','co2E_2')
basinArea_h<-pivot_wider(basinArea[,c('basinCode','Mon',cols)],
                         id_cols=basinCode,names_from=Mon,values_from=cols)
basinArea_h$basinCode<-as.character(basinArea_h$basinCode)
hydroBAS_res2<-left_join(hydroBAS,basinArea_h,by=c('basinCode'))
#correcting names
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('pco2_2')),funs(str_replace(.,'pco2_2','pco2')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('k_2')),funs(str_replace(.,'k_2','k')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('co2F_2')),funs(str_replace(.,'co2F_2','co2F')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('totArea_2')),funs(str_replace(.,'totArea_2','rivArea')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('totEphemArea_2')),funs(str_replace(.,'totEphemArea_2','ephemArea')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('totEffectArea_2')),funs(str_replace(.,'totEffectArea_2','effectArea')))
hydroBAS_res2<-hydroBAS_res2%>%rename_at(vars(starts_with('co2E_2')),funs(str_replace(.,'co2E_2','co2E')))
#using areaRatio to partition rivArea,ephemArea,icecovArea,and co2E
for(Mon in month.abb){hydroBAS_res2[,paste0('rivArea_',Mon)]<-hydroBAS_res2$areaRatio*hydroBAS_res2[,paste0('rivArea_',Mon)]}
for(Mon in month.abb){hydroBAS_res2[,paste0('ephemArea_',Mon)]<-hydroBAS_res2$areaRatio*hydroBAS_res2[,paste0('ephemArea_',Mon)]}
for(Mon in month.abb){hydroBAS_res2[,paste0('effectArea_',Mon)]<-hydroBAS_res2$areaRatio*hydroBAS_res2[,paste0('effectArea_',Mon)]}
for(Mon in month.abb){hydroBAS_res2[,paste0('co2E_',Mon)]<-hydroBAS_res2$areaRatio*hydroBAS_res2[,paste0('co2E_',Mon)]}
#correcting for ice covered area and co2E for part 2
hydroBAS_res2<-left_join(hydroBAS_res2,icecov,by='HYBAS_ID')
#calculate iceCovered area
hydroBAS_res2<-
  hydroBAS_res2%>%mutate(
    icecovArea_Jan=effectArea_Jan*iceCov_Jan,
    icecovArea_Feb=effectArea_Feb*iceCov_Feb,
    icecovArea_Mar=effectArea_Mar*iceCov_Mar,
    icecovArea_Apr=effectArea_Apr*iceCov_Apr,
    icecovArea_May=effectArea_May*iceCov_May,
    icecovArea_Jun=effectArea_Jun*iceCov_Jun,
    icecovArea_Jul=effectArea_Jul*iceCov_Jul,
    icecovArea_Aug=effectArea_Aug*iceCov_Aug,
    icecovArea_Sep=effectArea_Sep*iceCov_Sep,
    icecovArea_Oct=effectArea_Oct*iceCov_Oct,
    icecovArea_Nov=effectArea_Nov*iceCov_Nov,
    icecovArea_Dec=effectArea_Dec*iceCov_Dec,
    effectArea_Jan=effectArea_Jan-icecovArea_Jan,
    effectArea_Feb=effectArea_Feb-icecovArea_Feb,
    effectArea_Mar=effectArea_Mar-icecovArea_Mar,
    effectArea_Apr=effectArea_Apr-icecovArea_Apr,
    effectArea_May=effectArea_May-icecovArea_May,
    effectArea_Jun=effectArea_Jun-icecovArea_Jun,
    effectArea_Jul=effectArea_Jul-icecovArea_Jul,
    effectArea_Aug=effectArea_Aug-icecovArea_Aug,
    effectArea_Sep=effectArea_Sep-icecovArea_Sep,
    effectArea_Oct=effectArea_Oct-icecovArea_Oct,
    effectArea_Nov=effectArea_Nov-icecovArea_Nov,
    effectArea_Dec=effectArea_Dec-icecovArea_Dec,
    co2E_Jan=co2E_Jan*(1-iceCov_Jan),
    co2E_Feb=co2E_Feb*(1-iceCov_Feb),
    co2E_Mar=co2E_Mar*(1-iceCov_Mar),
    co2E_Apr=co2E_Apr*(1-iceCov_Apr),
    co2E_May=co2E_May*(1-iceCov_May),
    co2E_Jun=co2E_Jun*(1-iceCov_Jun),
    co2E_Jul=co2E_Jul*(1-iceCov_Jul),
    co2E_Aug=co2E_Aug*(1-iceCov_Aug),
    co2E_Sep=co2E_Sep*(1-iceCov_Sep),
    co2E_Oct=co2E_Oct*(1-iceCov_Oct),
    co2E_Nov=co2E_Nov*(1-iceCov_Nov),
    co2E_Dec=co2E_Dec*(1-iceCov_Dec)
  )

#combining res1 and res2
#for rivArea, ephemArea, icecovArea, effect area and co2E, just sum up
for(Mon in month.abb){hydroBAS[,paste0('rivArea_',Mon)]<-hydroBAS_res1[,paste0('rivArea_',Mon)]+hydroBAS_res2[,paste0('rivArea_',Mon)]}
for(Mon in month.abb){hydroBAS[,paste0('ephemArea_',Mon)]<-hydroBAS_res1[,paste0('ephemArea_',Mon)]+hydroBAS_res2[,paste0('ephemArea_',Mon)]}
for(Mon in month.abb){hydroBAS[,paste0('icecovArea_',Mon)]<-hydroBAS_res1[,paste0('icecovArea_',Mon)]+hydroBAS_res2[,paste0('icecovArea_',Mon)]}
for(Mon in month.abb){hydroBAS[,paste0('effectArea_',Mon)]<-hydroBAS_res1[,paste0('effectArea_',Mon)]+hydroBAS_res2[,paste0('effectArea_',Mon)]}
for(Mon in month.abb){hydroBAS[,paste0('co2E_',Mon)]<-hydroBAS_res1[,paste0('co2E_',Mon)]+hydroBAS_res2[,paste0('co2E_',Mon)]}

#for pco2,k and co2F, it weighted avgs
for(Mon in month.abb){
  hydroBAS[,paste0('pco2_',Mon)]<-(hydroBAS_res1[,paste0('effectArea_',Mon)]*hydroBAS_res1[,paste0('pco2_',Mon)]+
     hydroBAS_res2[,paste0('effectArea_',Mon)]*hydroBAS_res2[,paste0('pco2_',Mon)])/hydroBAS[,paste0('effectArea_',Mon)]
}
for(Mon in month.abb){hydroBAS[,paste0('k_',Mon)]<-
  (hydroBAS_res1[,paste0('effectArea_',Mon)]*hydroBAS_res1[,paste0('k_',Mon)]+
     hydroBAS_res2[,paste0('effectArea_',Mon)]*hydroBAS_res2[,paste0('k_',Mon)])/hydroBAS[,paste0('effectArea_',Mon)]}
for(Mon in month.abb){hydroBAS[,paste0('co2F_',Mon)]<-
  (hydroBAS_res1[,paste0('effectArea_',Mon)]*hydroBAS_res1[,paste0('co2F_',Mon)]+
     hydroBAS_res2[,paste0('effectArea_',Mon)]*hydroBAS_res2[,paste0('co2F_',Mon)])/hydroBAS[,paste0('effectArea_',Mon)]}
#correcting for basins where there is no effective area
for(Mon in month.abb){
  hydroBAS[hydroBAS[,paste0('effectArea_',Mon)]==0,paste0('pco2_',Mon)]<-0
  hydroBAS[hydroBAS[,paste0('effectArea_',Mon)]==0,paste0('k_',Mon)]<-0
  hydroBAS[hydroBAS[,paste0('effectArea_',Mon)]==0,paste0('co2F_',Mon)]<-0
  }
# rm(hydroBAS_res1,hydroBAS_res2,icecov,basinArea_h,GRWLwidth)
gc()

####stats by climate zones####
#1
hydroBASsum_1<-
  hydroBAS[!(hydroBAS$climzone%in%c('Polar')&hydroBAS$wetness%in%'Wet'),]%>%
  # group_by(climzone,wetness)%>%
  # group_by(climzone)%>%
  summarise(
    effectAreaSum_Jan=sum(effectArea_Jan),
    effectAreaSum_Feb=sum(effectArea_Feb),
    effectAreaSum_Mar=sum(effectArea_Mar),
    effectAreaSum_Apr=sum(effectArea_Apr),
    effectAreaSum_May=sum(effectArea_May),
    effectAreaSum_Jun=sum(effectArea_Jun),
    effectAreaSum_Jul=sum(effectArea_Jul),
    effectAreaSum_Aug=sum(effectArea_Aug),
    effectAreaSum_Sep=sum(effectArea_Sep),
    effectAreaSum_Oct=sum(effectArea_Oct),
    effectAreaSum_Nov=sum(effectArea_Nov),
    effectAreaSum_Dec=sum(effectArea_Dec),
    pco2_Jan=sum(pco2_Jan*effectArea_Jan)/effectAreaSum_Jan,
    pco2_Feb=sum(pco2_Feb*effectArea_Feb)/effectAreaSum_Feb,
    pco2_Mar=sum(pco2_Mar*effectArea_Mar)/effectAreaSum_Mar,
    pco2_Apr=sum(pco2_Apr*effectArea_Apr)/effectAreaSum_Apr,
    pco2_May=sum(pco2_May*effectArea_May)/effectAreaSum_May,
    pco2_Jun=sum(pco2_Jun*effectArea_Jun)/effectAreaSum_Jun,
    pco2_Jul=sum(pco2_Jul*effectArea_Jul)/effectAreaSum_Jul,
    pco2_Aug=sum(pco2_Aug*effectArea_Aug)/effectAreaSum_Aug,
    pco2_Sep=sum(pco2_Sep*effectArea_Sep)/effectAreaSum_Sep,
    pco2_Oct=sum(pco2_Oct*effectArea_Oct)/effectAreaSum_Oct,
    pco2_Nov=sum(pco2_Nov*effectArea_Nov)/effectAreaSum_Nov,
    pco2_Dec=sum(pco2_Dec*effectArea_Dec)/effectAreaSum_Dec,
    k_Jan=sum(k_Jan*effectArea_Jan)/effectAreaSum_Jan,
    k_Feb=sum(k_Feb*effectArea_Feb)/effectAreaSum_Feb,
    k_Mar=sum(k_Mar*effectArea_Mar)/effectAreaSum_Mar,
    k_Apr=sum(k_Apr*effectArea_Apr)/effectAreaSum_Apr,
    k_May=sum(k_May*effectArea_May)/effectAreaSum_May,
    k_Jun=sum(k_Jun*effectArea_Jun)/effectAreaSum_Jun,
    k_Jul=sum(k_Jul*effectArea_Jul)/effectAreaSum_Jul,
    k_Aug=sum(k_Aug*effectArea_Aug)/effectAreaSum_Aug,
    k_Sep=sum(k_Sep*effectArea_Sep)/effectAreaSum_Sep,
    k_Oct=sum(k_Oct*effectArea_Oct)/effectAreaSum_Oct,
    k_Nov=sum(k_Nov*effectArea_Nov)/effectAreaSum_Nov,
    k_Dec=sum(k_Dec*effectArea_Dec)/effectAreaSum_Dec,
    co2F_Jan=sum(co2F_Jan*effectArea_Jan)/effectAreaSum_Jan,
    co2F_Feb=sum(co2F_Feb*effectArea_Feb)/effectAreaSum_Feb,
    co2F_Mar=sum(co2F_Mar*effectArea_Mar)/effectAreaSum_Mar,
    co2F_Apr=sum(co2F_Apr*effectArea_Apr)/effectAreaSum_Apr,
    co2F_May=sum(co2F_May*effectArea_May)/effectAreaSum_May,
    co2F_Jun=sum(co2F_Jun*effectArea_Jun)/effectAreaSum_Jun,
    co2F_Jul=sum(co2F_Jul*effectArea_Jul)/effectAreaSum_Jul,
    co2F_Aug=sum(co2F_Aug*effectArea_Aug)/effectAreaSum_Aug,
    co2F_Sep=sum(co2F_Sep*effectArea_Sep)/effectAreaSum_Sep,
    co2F_Oct=sum(co2F_Oct*effectArea_Oct)/effectAreaSum_Oct,
    co2F_Nov=sum(co2F_Nov*effectArea_Nov)/effectAreaSum_Nov,
    co2F_Dec=sum(co2F_Dec*effectArea_Dec)/effectAreaSum_Dec)
# write_csv(hydroBASsum_1,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_1.csv'))
write_csv(hydroBASsum_1,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_global.csv'))


#2
hydroBASsum_2<-
  hydroBAS%>%
  # group_by(climzone,wetness)%>%
  group_by(climzone)%>%
  summarise(
    rivAreaSum_Jan=sum(rivArea_Jan),
    rivAreaSum_Feb=sum(rivArea_Feb),
    rivAreaSum_Mar=sum(rivArea_Mar),
    rivAreaSum_Apr=sum(rivArea_Apr),
    rivAreaSum_May=sum(rivArea_May),
    rivAreaSum_Jun=sum(rivArea_Jun),
    rivAreaSum_Jul=sum(rivArea_Jul),
    rivAreaSum_Aug=sum(rivArea_Aug),
    rivAreaSum_Sep=sum(rivArea_Sep),
    rivAreaSum_Oct=sum(rivArea_Oct),
    rivAreaSum_Nov=sum(rivArea_Nov),
    rivAreaSum_Dec=sum(rivArea_Dec),
    effectAreaSum_Jan=sum(effectArea_Jan),
    effectAreaSum_Feb=sum(effectArea_Feb),
    effectAreaSum_Mar=sum(effectArea_Mar),
    effectAreaSum_Apr=sum(effectArea_Apr),
    effectAreaSum_May=sum(effectArea_May),
    effectAreaSum_Jun=sum(effectArea_Jun),
    effectAreaSum_Jul=sum(effectArea_Jul),
    effectAreaSum_Aug=sum(effectArea_Aug),
    effectAreaSum_Sep=sum(effectArea_Sep),
    effectAreaSum_Oct=sum(effectArea_Oct),
    effectAreaSum_Nov=sum(effectArea_Nov),
    effectAreaSum_Dec=sum(effectArea_Dec),
    co2Esum_Jan=sum(co2E_Jan)/1000000,
    co2Esum_Feb=sum(co2E_Feb)/1000000,
    co2Esum_Mar=sum(co2E_Mar)/1000000,
    co2Esum_Apr=sum(co2E_Apr)/1000000,
    co2Esum_May=sum(co2E_May)/1000000,
    co2Esum_Jun=sum(co2E_Jun)/1000000,
    co2Esum_Jul=sum(co2E_Jul)/1000000,
    co2Esum_Aug=sum(co2E_Aug)/1000000,
    co2Esum_Sep=sum(co2E_Sep)/1000000,
    co2Esum_Oct=sum(co2E_Oct)/1000000,
    co2Esum_Nov=sum(co2E_Nov)/1000000,
    co2Esum_Dec=sum(co2E_Dec)/1000000,
    ephemAreaSum_Jan=sum(ephemArea_Jan),
    ephemAreaSum_Feb=sum(ephemArea_Feb),
    ephemAreaSum_Mar=sum(ephemArea_Mar),
    ephemAreaSum_Apr=sum(ephemArea_Apr),
    ephemAreaSum_May=sum(ephemArea_May),
    ephemAreaSum_Jun=sum(ephemArea_Jun),
    ephemAreaSum_Jul=sum(ephemArea_Jul),
    ephemAreaSum_Aug=sum(ephemArea_Aug),
    ephemAreaSum_Sep=sum(ephemArea_Sep),
    ephemAreaSum_Oct=sum(ephemArea_Oct),
    ephemAreaSum_Nov=sum(ephemArea_Nov),
    ephemAreaSum_Dec=sum(ephemArea_Dec),
    icecovAreaSum_Jan=sum(icecovArea_Jan),
    icecovAreaSum_Feb=sum(icecovArea_Feb),
    icecovAreaSum_Mar=sum(icecovArea_Mar),
    icecovAreaSum_Apr=sum(icecovArea_Apr),
    icecovAreaSum_May=sum(icecovArea_May),
    icecovAreaSum_Jun=sum(icecovArea_Jun),
    icecovAreaSum_Jul=sum(icecovArea_Jul),
    icecovAreaSum_Aug=sum(icecovArea_Aug),
    icecovAreaSum_Sep=sum(icecovArea_Sep),
    icecovAreaSum_Oct=sum(icecovArea_Oct),
    icecovAreaSum_Nov=sum(icecovArea_Nov),
    icecovAreaSum_Dec=sum(icecovArea_Dec))

# write_csv(hydroBASsum_2,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_2.csv'))


#Amazon basins:
#South Amazon,6040262110; North Amazaon,6040280400
#Congo basins:1041213630,1041213640,1041174950,1041156960,1040020570
#up Congo, 1041156960; north congo,1041156950; south congo,1041213640.
#SE basins:5040000010,5040004220,5040015820,5040015830,5040021180,5040041880,5040046910,4040018350

hydroBASsum_s<-
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('6040262220','6040262110','6040285990','6040280410','6040280400'),]%>%
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('6040262110'),]%>% #south Amazon
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('6040280400'),]%>% #north Amazon
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('1041213630','1041213640','1041174950','1041156960','1040020570'),]%>%
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('1041156960'),]%>% #up Congo
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('1041156950'),]%>% #north Congo
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('1041213640'),]%>% #south Congo
  # hydroBAS[hydroBAS$HYBAS_ID%in%c('5040000010','5040004220','5040015820',
                                  # '5040015830','5040021180','5040041880','5040046910','4040018350'),]%>%
  summarise(
    effectAreaSum_Jan=sum(effectArea_Jan),
    effectAreaSum_Feb=sum(effectArea_Feb),
    effectAreaSum_Mar=sum(effectArea_Mar),
    effectAreaSum_Apr=sum(effectArea_Apr),
    effectAreaSum_May=sum(effectArea_May),
    effectAreaSum_Jun=sum(effectArea_Jun),
    effectAreaSum_Jul=sum(effectArea_Jul),
    effectAreaSum_Aug=sum(effectArea_Aug),
    effectAreaSum_Sep=sum(effectArea_Sep),
    effectAreaSum_Oct=sum(effectArea_Oct),
    effectAreaSum_Nov=sum(effectArea_Nov),
    effectAreaSum_Dec=sum(effectArea_Dec),
    pco2_Jan=sum(pco2_Jan*effectArea_Jan)/effectAreaSum_Jan,
    pco2_Feb=sum(pco2_Feb*effectArea_Feb)/effectAreaSum_Feb,
    pco2_Mar=sum(pco2_Mar*effectArea_Mar)/effectAreaSum_Mar,
    pco2_Apr=sum(pco2_Apr*effectArea_Apr)/effectAreaSum_Apr,
    pco2_May=sum(pco2_May*effectArea_May)/effectAreaSum_May,
    pco2_Jun=sum(pco2_Jun*effectArea_Jun)/effectAreaSum_Jun,
    pco2_Jul=sum(pco2_Jul*effectArea_Jul)/effectAreaSum_Jul,
    pco2_Aug=sum(pco2_Aug*effectArea_Aug)/effectAreaSum_Aug,
    pco2_Sep=sum(pco2_Sep*effectArea_Sep)/effectAreaSum_Sep,
    pco2_Oct=sum(pco2_Oct*effectArea_Oct)/effectAreaSum_Oct,
    pco2_Nov=sum(pco2_Nov*effectArea_Nov)/effectAreaSum_Nov,
    pco2_Dec=sum(pco2_Dec*effectArea_Dec)/effectAreaSum_Dec,
    k_Jan=sum(k_Jan*effectArea_Jan)/effectAreaSum_Jan,
    k_Feb=sum(k_Feb*effectArea_Feb)/effectAreaSum_Feb,
    k_Mar=sum(k_Mar*effectArea_Mar)/effectAreaSum_Mar,
    k_Apr=sum(k_Apr*effectArea_Apr)/effectAreaSum_Apr,
    k_May=sum(k_May*effectArea_May)/effectAreaSum_May,
    k_Jun=sum(k_Jun*effectArea_Jun)/effectAreaSum_Jun,
    k_Jul=sum(k_Jul*effectArea_Jul)/effectAreaSum_Jul,
    k_Aug=sum(k_Aug*effectArea_Aug)/effectAreaSum_Aug,
    k_Sep=sum(k_Sep*effectArea_Sep)/effectAreaSum_Sep,
    k_Oct=sum(k_Oct*effectArea_Oct)/effectAreaSum_Oct,
    k_Nov=sum(k_Nov*effectArea_Nov)/effectAreaSum_Nov,
    k_Dec=sum(k_Dec*effectArea_Dec)/effectAreaSum_Dec,
    co2F_Jan=sum(co2F_Jan*effectArea_Jan)/effectAreaSum_Jan,
    co2F_Feb=sum(co2F_Feb*effectArea_Feb)/effectAreaSum_Feb,
    co2F_Mar=sum(co2F_Mar*effectArea_Mar)/effectAreaSum_Mar,
    co2F_Apr=sum(co2F_Apr*effectArea_Apr)/effectAreaSum_Apr,
    co2F_May=sum(co2F_May*effectArea_May)/effectAreaSum_May,
    co2F_Jun=sum(co2F_Jun*effectArea_Jun)/effectAreaSum_Jun,
    co2F_Jul=sum(co2F_Jul*effectArea_Jul)/effectAreaSum_Jul,
    co2F_Aug=sum(co2F_Aug*effectArea_Aug)/effectAreaSum_Aug,
    co2F_Sep=sum(co2F_Sep*effectArea_Sep)/effectAreaSum_Sep,
    co2F_Oct=sum(co2F_Oct*effectArea_Oct)/effectAreaSum_Oct,
    co2F_Nov=sum(co2F_Nov*effectArea_Nov)/effectAreaSum_Nov,
    co2F_Dec=sum(co2F_Dec*effectArea_Dec)/effectAreaSum_Dec,
    rivAreaSum_Jan=sum(rivArea_Jan),
    rivAreaSum_Feb=sum(rivArea_Feb),
    rivAreaSum_Mar=sum(rivArea_Mar),
    rivAreaSum_Apr=sum(rivArea_Apr),
    rivAreaSum_May=sum(rivArea_May),
    rivAreaSum_Jun=sum(rivArea_Jun),
    rivAreaSum_Jul=sum(rivArea_Jul),
    rivAreaSum_Aug=sum(rivArea_Aug),
    rivAreaSum_Sep=sum(rivArea_Sep),
    rivAreaSum_Oct=sum(rivArea_Oct),
    rivAreaSum_Nov=sum(rivArea_Nov),
    rivAreaSum_Dec=sum(rivArea_Dec),
    effectAreaSum_Jan=sum(effectArea_Jan),
    effectAreaSum_Feb=sum(effectArea_Feb),
    effectAreaSum_Mar=sum(effectArea_Mar),
    effectAreaSum_Apr=sum(effectArea_Apr),
    effectAreaSum_May=sum(effectArea_May),
    effectAreaSum_Jun=sum(effectArea_Jun),
    effectAreaSum_Jul=sum(effectArea_Jul),
    effectAreaSum_Aug=sum(effectArea_Aug),
    effectAreaSum_Sep=sum(effectArea_Sep),
    effectAreaSum_Oct=sum(effectArea_Oct),
    effectAreaSum_Nov=sum(effectArea_Nov),
    effectAreaSum_Dec=sum(effectArea_Dec),
    co2Esum_Jan=sum(co2E_Jan)/1000000,
    co2Esum_Feb=sum(co2E_Feb)/1000000,
    co2Esum_Mar=sum(co2E_Mar)/1000000,
    co2Esum_Apr=sum(co2E_Apr)/1000000,
    co2Esum_May=sum(co2E_May)/1000000,
    co2Esum_Jun=sum(co2E_Jun)/1000000,
    co2Esum_Jul=sum(co2E_Jul)/1000000,
    co2Esum_Aug=sum(co2E_Aug)/1000000,
    co2Esum_Sep=sum(co2E_Sep)/1000000,
    co2Esum_Oct=sum(co2E_Oct)/1000000,
    co2Esum_Nov=sum(co2E_Nov)/1000000,
    co2Esum_Dec=sum(co2E_Dec)/1000000,
    ephemAreaSum_Jan=sum(ephemArea_Jan),
    ephemAreaSum_Feb=sum(ephemArea_Feb),
    ephemAreaSum_Mar=sum(ephemArea_Mar),
    ephemAreaSum_Apr=sum(ephemArea_Apr),
    ephemAreaSum_May=sum(ephemArea_May),
    ephemAreaSum_Jun=sum(ephemArea_Jun),
    ephemAreaSum_Jul=sum(ephemArea_Jul),
    ephemAreaSum_Aug=sum(ephemArea_Aug),
    ephemAreaSum_Sep=sum(ephemArea_Sep),
    ephemAreaSum_Oct=sum(ephemArea_Oct),
    ephemAreaSum_Nov=sum(ephemArea_Nov),
    ephemAreaSum_Dec=sum(ephemArea_Dec),
    icecovAreaSum_Jan=sum(icecovArea_Jan),
    icecovAreaSum_Feb=sum(icecovArea_Feb),
    icecovAreaSum_Mar=sum(icecovArea_Mar),
    icecovAreaSum_Apr=sum(icecovArea_Apr),
    icecovAreaSum_May=sum(icecovArea_May),
    icecovAreaSum_Jun=sum(icecovArea_Jun),
    icecovAreaSum_Jul=sum(icecovArea_Jul),
    icecovAreaSum_Aug=sum(icecovArea_Aug),
    icecovAreaSum_Sep=sum(icecovArea_Sep),
    icecovAreaSum_Oct=sum(icecovArea_Oct),
    icecovAreaSum_Nov=sum(icecovArea_Nov),
    icecovAreaSum_Dec=sum(icecovArea_Dec))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_Amazon.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_sAmazon.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_nAmazon.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_Africa.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_uAfrica.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_nAfrica.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_sAfrica.csv'))
# write_csv(hydroBASsum_s,paste0(wd,'/output/table/regionSurfArea/hydroBASsum_SE.csv'))


# #plotting
# hydroBASsum_1[,c('climzone','wetness',paste0('pco2_',month.abb))]%>%
#   pivot_longer(cols=c(-climzone,-wetness),names_to='Mon',values_to='pco2')%>%
#   separate(Mon,into=c('tp','Mon'))%>%
#   ggplot(aes(x=factor(Mon,levels=month.abb),y=pco2,group=1))+
#   geom_point(size=2,color='red')+
#   geom_line()+
#   labs(x='Mon',y=expression(italic(p)*'C'*O[2]*' ('*mu*'atm)'))+
#   scale_x_discrete(breaks=c('Feb','Apr','Jun','Aug','Oct','Dec'))+
#   scale_y_continuous(limits=c(300,3800),breaks=seq(1000,3000,by=1000))+
#   theme_classic()+
#   theme(panel.border=element_rect(fill=NA,size=0.5),
#         strip.background=element_blank())+
#   facet_grid(climzone~wetness)
# ggsave(paste0(wd,'/output/figure/mainfigs/pco2_climzone.png'),
#        width=20,height=18,units='cm')

# hydroBASsum_1[,c('climzone','wetness',paste0('k_',month.abb))]%>%
#   pivot_longer(cols=c(-climzone,-wetness),names_to='Mon',values_to='k')%>%
#   separate(Mon,into=c('tp','Mon'))%>%
#   ggplot(aes(x=factor(Mon,levels=month.abb),y=k,group=1))+
#   geom_point(size=2,color='red')+
#   geom_line()+
#   labs(x='Mon',y=expression('k (m '*d^-1*')'))+
#   scale_x_discrete(breaks=c('Feb','Apr','Jun','Aug','Oct','Dec'))+
#   scale_y_continuous(limits=c(0,35),breaks=seq(0,30,by=10))+
#   theme_classic()+
#   theme(panel.border=element_rect(fill=NA,size=0.5),
#         strip.background=element_blank())+
#   facet_grid(climzone~wetness)
# ggsave(paste0(wd,'/output/figure/mainfigs/k_climzone.png'),
#        width=20,height=18,units='cm')
# 
# hydroBASsum_1[,c('climzone','wetness',paste0('co2F_',month.abb))]%>%
#   pivot_longer(cols=c(-climzone,-wetness),names_to='Mon',values_to='co2F')%>%
#   separate(Mon,into=c('tp','Mon'))%>%
#   ggplot(aes(x=factor(Mon,levels=month.abb),y=co2F,group=1))+
#   geom_point(size=2,color='red')+
#   geom_line()+
#   labs(x='Mon',y=expression('C'*O[2]*' Efflux (g C '*m^-2*' '*yr^-1*')'))+
#   scale_x_discrete(breaks=c('Feb','Apr','Jun','Aug','Oct','Dec'))+
#   scale_y_continuous(limits=c(0,7000),breaks=seq(1000,5000,by=2000))+
#   theme_classic()+
#   theme(panel.border=element_rect(fill=NA,size=0.5),
#         strip.background=element_blank())+
#   facet_grid(climzone~wetness)
# ggsave(paste0(wd,'/output/figure/mainfigs/co2F_climzone.png'),
#        width=20,height=18,units='cm')

####comparing to major terrestrial fluxes####
#check units of GPP, NPP and SR
sum(hydroBAS$GPP*hydroBAS$basinArea/1000000000,na.rm=TRUE) #118 PgYr-1
sum(hydroBAS$NPP*hydroBAS$basinArea/1000000000,na.rm=TRUE) #56.9 PgYr-1
sum(hydroBAS$pyearRS*hydroBAS$basinArea/1000000000,na.rm=TRUE) #96 PgYr-1
sum(hydroBAS[,paste0('npp_',str_pad(1:12,2,pad='0'))]%>%rowMeans()
    *hydroBAS$basinArea/1000000000,na.rm=TRUE)#73.23
sum(hydroBAS[,paste0('gpp_',str_pad(1:12,2,pad='0'))]%>%rowMeans()
    *hydroBAS$basinArea/1000000000,na.rm=TRUE)#119.67

sum(hydroBAS$pRS_12*hydroBAS$basinArea/1000000000,na.rm=TRUE)
sum(hydroBAS$gpp_12*hydroBAS$basinArea/1000000000,na.rm=TRUE)
sum(hydroBAS$npp_01*hydroBAS$basinArea/1000000000,na.rm=TRUE)

#calc %rivArea
hydroBAS[,paste0('effectAr_Bas_',month.abb)]<-
  hydroBAS[,paste0('effectArea_',month.abb)]/hydroBAS$basinArea
hydroBAS$effectArea_Ann<-hydroBAS[,paste0('effectArea_',month.abb)]%>%rowMeans()
hydroBAS$effectAr_Bas_Ann<-hydroBAS$effectArea_Ann/hydroBAS$basinArea
#some stats
sum(hydroBAS$effectArea_Ann)/sum(hydroBAS$basinArea)*100#~0.50% of the land surface area is streams and rivers
sum(hydroBAS[hydroBAS$Lat>=0,]$basinArea)/sum(hydroBAS$basinArea)*100 #73.6% of land surface is northern hemisphere

#calc GPP, NPP and SR at the basin scale
hydroBAS[,paste0('gpp_Bas_',str_pad(1:12,2,pad='0'))]<-
  hydroBAS[,paste0('gpp_',str_pad(1:12,2,pad='0'))]*hydroBAS$basinArea/1000 #10^9 gCyr-1
hydroBAS[,paste0('npp_Bas_',str_pad(1:12,2,pad='0'))]<-
  abs(hydroBAS[,paste0('npp_',str_pad(1:12,2,pad='0'))])*hydroBAS$basinArea/1000 #10^9 gCyr-1
hydroBAS[,paste0('SR_Bas_',str_pad(1:12,2,pad='0'))]<-
  hydroBAS[,paste0('pRS_',str_pad(1:12,2,pad='0'))]*hydroBAS$basinArea/1000 #10^9 gCyr-1

#calc ratios of co2E to GPP, NPP, and SR
#co2r GPP
hydroBAS[,paste0('co2r_gpp_',month.abb)]<-
  hydroBAS[,paste0('co2E_',month.abb)]/hydroBAS[,paste0('gpp_Bas_',str_pad(1:12,2,pad='0'))]
hydroBAS$co2r_gpp_Ann<-hydroBAS[,paste0('co2r_gpp_',month.abb)]%>%rowMeans()
#co2r NPP
hydroBAS[,paste0('co2r_npp_',month.abb)]<-
  hydroBAS[,paste0('co2E_',month.abb)]/hydroBAS[,paste0('npp_Bas_',str_pad(1:12,2,pad='0'))]
hydroBAS$co2r_npp_Ann<-hydroBAS[,paste0('co2r_npp_',month.abb)]%>%rowMeans()
#co2r SR
hydroBAS[,paste0('co2r_SR_',month.abb)]<-
  hydroBAS[,paste0('co2E_',month.abb)]/hydroBAS[,paste0('SR_Bas_',str_pad(1:12,2,pad='0'))]
hydroBAS$co2r_SR_Ann<-hydroBAS[,paste0('co2r_SR_',month.abb)]%>%rowMeans()
#exporting results for mapping
write_csv(hydroBAS[hydroBAS$runoff>0&hydroBAS$basinArea>2000,
                   c('HYBAS_ID','effectAr_Bas_Ann','co2r_gpp_Ann','co2r_npp_Ann','co2r_SR_Ann','runoff')],
          paste0(wd,'/output/table/regionSurfArea/co2Eratios.csv'))

co2Eratios<-read_csv(paste0(wd,'/output/table/regionSurfArea/co2Eratios.csv'),
                     col_types=cols(.default='d',HYBAS_ID='c'))
sum(is.na(co2Eratios$co2r_gpp_Ann))
sum(is.infinite(co2Eratios$co2r_gpp_Ann))
quantile(co2Eratios$co2r_gpp_Ann,0.8,na.rm=TRUE)
co2Eratios[is.infinite(co2Eratios$co2r_gpp_Ann),]$co2r_gpp_Ann<-0.8
sum(is.na(co2Eratios$co2r_npp_Ann))
sum(is.infinite(co2Eratios$co2r_npp_Ann))
# co2Eratios[is.na(co2Eratios$co2r_SR_Ann),]$co2r_SR_Ann<-0
# co2Eratios[is.infinite(co2Eratios$co2r_SR_Ann),]$co2r_SR_Ann<-0
sum(is.na(co2Eratios$co2r_SR_Ann))
sum(is.infinite(co2Eratios$co2r_SR_Ann))
names(co2Eratios)<-c('HYBAS_ID','Ar','co2rgpp','co2rnpp','co2rSR','runoff')
# write_csv(co2Eratios,paste0(wd,'/output/table/regionSurfArea/co2Eratios.csv'))
co2Eratios[co2Eratios$HYBAS_ID=='1040040510','co2rSR']<-0.008
co2EratiosSR<-co2Eratios[!is.na(co2Eratios$co2rSR),]
write_csv(co2EratiosSR[co2EratiosSR$runoff>=0,],
          paste0(wd,'/output/table/regionSurfArea/co2Eratios_sr.csv'))

#stats by climate zones
hydroBAS_ratios<-
  # hydroBAS[!is.na(hydroBAS$pRS_01),]%>%group_by(climzone)%>%summarise(
  hydroBAS[hydroBAS$runoff>15&hydroBAS$basinArea>2000,]%>%group_by(wetness)%>%summarise(
  # hydroBAS[hydroBAS$runoff>15&hydroBAS$basinArea>2000,]%>%group_by(climzone,wetness)%>%summarise(
    effAr_Bas_Jan=sum(effectArea_Jan)/sum(basinArea)*100,
    effAr_Bas_Feb=sum(effectArea_Feb)/sum(basinArea)*100,
    effAr_Bas_Mar=sum(effectArea_Mar)/sum(basinArea)*100,
    effAr_Bas_Apr=sum(effectArea_Apr)/sum(basinArea)*100,
    effAr_Bas_May=sum(effectArea_May)/sum(basinArea)*100,
    effAr_Bas_Jun=sum(effectArea_Jun)/sum(basinArea)*100,
    effAr_Bas_Jul=sum(effectArea_Jul)/sum(basinArea)*100,
    effAr_Bas_Aug=sum(effectArea_Aug)/sum(basinArea)*100,
    effAr_Bas_Sep=sum(effectArea_Sep)/sum(basinArea)*100,
    effAr_Bas_Oct=sum(effectArea_Oct)/sum(basinArea)*100,
    effAr_Bas_Nov=sum(effectArea_Nov)/sum(basinArea)*100,
    effAr_Bas_Dec=sum(effectArea_Dec)/sum(basinArea)*100,
    co2r_gpp_Jan=sum(co2E_Jan)/sum(gpp_Bas_01,na.rm=TRUE)*100,
    co2r_gpp_Feb=sum(co2E_Feb)/sum(gpp_Bas_02,na.rm=TRUE)*100,
    co2r_gpp_Mar=sum(co2E_Mar)/sum(gpp_Bas_03,na.rm=TRUE)*100,
    co2r_gpp_Apr=sum(co2E_Apr)/sum(gpp_Bas_04,na.rm=TRUE)*100,
    co2r_gpp_May=sum(co2E_May)/sum(gpp_Bas_05,na.rm=TRUE)*100,
    co2r_gpp_Jun=sum(co2E_Jun)/sum(gpp_Bas_06,na.rm=TRUE)*100,
    co2r_gpp_Jul=sum(co2E_Jul)/sum(gpp_Bas_07,na.rm=TRUE)*100,
    co2r_gpp_Aug=sum(co2E_Aug)/sum(gpp_Bas_08,na.rm=TRUE)*100,
    co2r_gpp_Sep=sum(co2E_Sep)/sum(gpp_Bas_09,na.rm=TRUE)*100,
    co2r_gpp_Oct=sum(co2E_Oct)/sum(gpp_Bas_10,na.rm=TRUE)*100,
    co2r_gpp_Nov=sum(co2E_Nov)/sum(gpp_Bas_11,na.rm=TRUE)*100,
    co2r_gpp_Dec=sum(co2E_Dec)/sum(gpp_Bas_12,na.rm=TRUE)*100,
    co2r_npp_Jan=sum(co2E_Jan)/sum(npp_Bas_01,na.rm=TRUE)*100,
    co2r_npp_Feb=sum(co2E_Feb)/sum(npp_Bas_02,na.rm=TRUE)*100,
    co2r_npp_Mar=sum(co2E_Mar)/sum(npp_Bas_03,na.rm=TRUE)*100,
    co2r_npp_Apr=sum(co2E_Apr)/sum(npp_Bas_04,na.rm=TRUE)*100,
    co2r_npp_May=sum(co2E_May)/sum(npp_Bas_05,na.rm=TRUE)*100,
    co2r_npp_Jun=sum(co2E_Jun)/sum(npp_Bas_06,na.rm=TRUE)*100,
    co2r_npp_Jul=sum(co2E_Jul)/sum(npp_Bas_07,na.rm=TRUE)*100,
    co2r_npp_Aug=sum(co2E_Aug)/sum(npp_Bas_08,na.rm=TRUE)*100,
    co2r_npp_Sep=sum(co2E_Sep)/sum(npp_Bas_09,na.rm=TRUE)*100,
    co2r_npp_Oct=sum(co2E_Oct)/sum(npp_Bas_10,na.rm=TRUE)*100,
    co2r_npp_Nov=sum(co2E_Nov)/sum(npp_Bas_11,na.rm=TRUE)*100,
    co2r_npp_Dec=sum(co2E_Dec)/sum(npp_Bas_12,na.rm=TRUE)*100,
    co2r_SR_Jan=sum(co2E_Jan)/sum(SR_Bas_01,na.rm=TRUE)*100,
    co2r_SR_Feb=sum(co2E_Feb)/sum(SR_Bas_02,na.rm=TRUE)*100,
    co2r_SR_Mar=sum(co2E_Mar)/sum(SR_Bas_03,na.rm=TRUE)*100,
    co2r_SR_Apr=sum(co2E_Apr)/sum(SR_Bas_04,na.rm=TRUE)*100,
    co2r_SR_May=sum(co2E_May)/sum(SR_Bas_05,na.rm=TRUE)*100,
    co2r_SR_Jun=sum(co2E_Jun)/sum(SR_Bas_06,na.rm=TRUE)*100,
    co2r_SR_Jul=sum(co2E_Jul)/sum(SR_Bas_07,na.rm=TRUE)*100,
    co2r_SR_Aug=sum(co2E_Aug)/sum(SR_Bas_08,na.rm=TRUE)*100,
    co2r_SR_Sep=sum(co2E_Sep)/sum(SR_Bas_09,na.rm=TRUE)*100,
    co2r_SR_Oct=sum(co2E_Oct)/sum(SR_Bas_10,na.rm=TRUE)*100,
    co2r_SR_Nov=sum(co2E_Nov)/sum(SR_Bas_11,na.rm=TRUE)*100,
    co2r_SR_Dec=sum(co2E_Dec)/sum(SR_Bas_12,na.rm=TRUE)*100
)

hydroBAS_ratios$wetness<-c('Arid','Moderate','Wet')

# pivot the dataset
hydroBAS_ratios<-
  # hydroBAS_ratios%>%pivot_longer(cols=c(-climzone),names_to='Var',values_to='Ratio')
  hydroBAS_ratios%>%pivot_longer(cols=c(-wetness),names_to='Var',values_to='Ratio')
  # hydroBAS_ratios%>%pivot_longer(cols=c(-climzone,-wetness),names_to='Var',values_to='Ratio')
hydroBAS_ratios<-
  hydroBAS_ratios%>%separate(Var,into=c('Var1','Var2','Mon'))
hydroBAS_ratios$Mon<-factor(hydroBAS_ratios$Mon,levels=month.abb)
hydroBAS_ratios$Var2<-
  factor(hydroBAS_ratios$Var2,levels=c('Bas','gpp','npp','SR'),
         labels=c('% River Surface Area','CO2 Emission : GPP',
                  'CO2 Emission : NPP', 'CO2 Emission : SR'))
ratioSum<-
  hydroBAS_ratios%>%group_by(wetness,Var2)%>%summarise(
  mn=mean(Ratio),
  rgmn=min(Ratio),
  rgmx=max(Ratio))#%>%
  # pivot_wider(id_cols=wetness,names_from=Var2,values_from=c('mn','rgmn','rgmx'))
write_csv(ratioSum,paste0(wd,'/output/table/regionSurfArea/ratioSum.csv'))

#calc annMean of co2E
hydroBAS$co2E_Ann<-hydroBAS[,paste0('co2E_',month.abb)]%>%rowMeans()
#normalizing co2E to basinArea
hydroBAS$co2E_Ann_Bas<-hydroBAS$co2E_Ann/hydroBAS$basinArea*1000 #gCm-2yr-1
for(Mon in month.abb){hydroBAS[,paste0('co2E_',Mon,'_Bas')]<-hydroBAS[,paste0('co2E_',Mon)]/hydroBAS$basinArea*1000}

#% River Surface Area versus carbon emission
hydroBAS_ratios_v<-
  hydroBAS[hydroBAS$basinArea>2000&hydroBAS$runoff>15&hydroBAS$GPP>500,c('climzone','runoff','effectAr_Bas_Ann','co2r_gpp_Ann','co2r_npp_Ann','co2r_SR_Ann')]%>%
  gather(key='Var',value='Ratio',-effectAr_Bas_Ann,-climzone,-runoff)
hydroBAS_ratios_v<-hydroBAS_ratios_v[!is.na(hydroBAS_ratios_v$Ratio),]
hydroBAS_ratios_v$Var<-
  factor(hydroBAS_ratios_v$Var,levels=c('co2r_gpp_Ann','co2r_npp_Ann','co2r_SR_Ann'),
         labels=c('CO2 Emission : GPP','CO2 Emission : NPP','CO2 Emission : Soil Resp.'))

hydroBAS_ratios_v<-
  hydroBAS_ratios_v[hydroBAS_ratios_v$Var%in%c('CO2 Emission : Soil Resp.'),]

#2
pB<-
  hydroBAS_ratios[hydroBAS_ratios$Var2%in%c('CO2 Emission : SR'),]%>%
  ggplot(aes(Mon,Ratio,color=wetness,group=wetness))+
  geom_point(size=3,fill='black',shape=1,stroke=1.2)+
  geom_line(size=0.6)+
  labs(x=expression('Mon'),
       y='% Soil Resp.')+
  scale_x_discrete(breaks=c('Feb','Apr','Jun','Aug','Oct','Dec'))+
  scale_y_continuous(limits=c(0,6),
                     breaks=seq(0,6,2))+
  # scale_y_continuous(limits=c(0.03,5), trans=log_trans(base=10),
  #                    breaks=trans_breaks('log10',function(x) 10^x)(c(0.1,5)),
  #                    labels=function(x){format(x,digit=1,drop0trailing=TRUE)})+
  scale_color_viridis(discrete=TRUE)+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.5,0.9),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.width=unit(2,'line'),
        legend.key.height=unit(1.2,'line'),
        legend.direction='horizontal',
        axis.ticks.length=unit(0.2,'cm'),
        strip.background=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
  )

pC<-
  hydroBAS_ratios_v[hydroBAS_ratios_v$Ratio<1,]%>%
  ggplot(aes(effectAr_Bas_Ann*100,Ratio*100,color=climzone))+
  geom_point(size=2.5,alpha=0.5)+
  scale_x_continuous(
    trans=log_trans(base=10),
    breaks=trans_breaks('log10',function(x) 10^x,n=3),
    labels=function(x){format(x,scientific=FALSE,drop0trailing=T)}
    )+
  scale_y_continuous(
    trans=log_trans(base=10),
    breaks=trans_breaks('log10',function(x) 10^x,n=3),
    labels=function(x){format(x,scientific=FALSE,drop0trailing=T)})+
  labs(x=expression('% River Surface Area'),
       y=expression('% Soil Resp.'),
       color=expression('Climatic Zone'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        strip.background=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.position=c(0.72,0.22))

ggsave(paste0(wd,'/output/figure/mainfigs/co2r.png'),
       plot=ggarrange(pB,pC,nrow=1,labels=c('',''),heights=c(1,1)),
       width=20,height=10,units='cm')
