library(foreign)
library(ncdf4)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(scales)
library(readxl)

wd <- 'D:/Research/globalEmission'

#Common variables, function, and datasets
Qquantiles <- c('Q0', 'Q10', 'Q20', 'Q30', 'Q40', 'Q50', 
                'Q60', 'Q70', 'Q80', 'Q90', 'Q100', 'Qmean', 'Qstd')

#monthly q stats extractor
qStatsExt <- function(nc,statsNm){
  qStats <- ncvar_get(nc, statsNm)
  qStats[qStats < 0.000000001]=0.000000001
  qStats <- data.frame(qStats)
  colnames(qStats) <- paste0(month.abb,statsNm) 
  if (nrow(qStats) == nrow(df)){
    return(qStats)
  } else {print('number of rows do not match ...')}
}
#ann q stat extractor
qAnnStatsExt<-function(nc,statsNm){
  qStats<-ncvar_get(nc,statsNm)
  qStats[qStats < 0.000000001]=0.000000001
  qStats<-data.frame(qStats)
  colnames(qStats)<-statsNm
  if (nrow(qStats) == nrow(df)){
    return(qStats)
  } else {print('number of rows do not match ...')}
}
#air-water temperature transformer
airwaterT<-function(ta){tw=0.67*ta+7.45}

for(i in 1:8){
  print(paste0('Processing ',i))
  df=read.dbf(paste0(wd,'/dataset/MeritHydroLinPan/pfaf_0',i,'_riv_3sMERIT.dbf'))
  df<-df[,c('COMID','strmOrder','Slope','Length')]
  
  qAnnStatsNC<-nc_open(paste0(wd,'/dataset/qStat/pfaf_0',i,'_all.nc'))
  df<-cbind(df,qAnnStatsExt(qAnnStatsNC,Qquantiles[12]))
  names(df)[5]<-'yeaQmean'
  
  #join monthly Q
  qStatsNC <- nc_open(paste0(wd, '/dataset/qStat/pfaf_0',i,'.nc'))
  df<-cbind(df,qStatsExt(qStatsNC, Qquantiles[12]))
  nc_close(qStatsNC)
  
  #calc month V and eD
  #averaging V as in Raymond2013
  # df[,c(paste0(month.abb,'V'))]<-data.frame(sapply(df[,c(paste0(month.abb,'Qmean'))],
  #                                                  function(x){exp(0.12*log(x)-1.06)/2+exp(0.285*log(x)-1.64)/2}))
  df[,c(paste0(month.abb,'V'))]<-data.frame(sapply(df[,c(paste0(month.abb,'Qmean'))],
                                                   function(x){exp(0.12*log(x)-1.06)}))#using USGS Q-V relationship
  df[,paste0(month.abb,'_eD')]<-data.frame(sapply(df[,paste0(month.abb,'V')],function(x){9.81*df$Slope*x}))

  #join temperature
  temp<-read_csv(paste0('D:/Research/globalEmission/output/table/MeritHydro/monTemp_0',i,'.csv'))
  colnames(temp)[2:13]<-paste0(month.abb,'_Ta')
  #calc Tw, correct Tw<0, calc Schmidt no and join Tw, Schmidt
  temp[,paste0(month.abb,'_Tw')]<-sapply(temp[2:13],airwaterT)
  for(col in paste0(month.abb,'_Tw')){temp[temp[,col]<0,col]=0}
  temp[,paste0(month.abb,'_Sc')]<-sapply(temp[,paste0(month.abb,'_Tw')],function(x){1742-91.24*x+2.208*(x^2)-0.0219*(x^3)})
  df<-left_join(df,temp[,c('COMID',paste0(month.abb,'_Tw'),paste0(month.abb,'_Sc'))],by='COMID')
  
  #join elev
  elev<-read_csv(paste0(wd,'/output/table/MeritHydro/eleSlope_0',i,'.csv'))
  elev<-elev[,c('COMID','elev')]
  df<-left_join(df,elev,by='COMID')

  df_1<-df[df$Slope<=0.01,]
  df_1<-df_1%>%mutate(
    Jan_k=(2841*Slope*JanV+2.02)*(600/Jan_Sc)^0.5,
    Feb_k=(2841*Slope*FebV+2.02)*(600/Feb_Sc)^0.5,
    Mar_k=(2841*Slope*MarV+2.02)*(600/Mar_Sc)^0.5,
    Apr_k=(2841*Slope*AprV+2.02)*(600/Apr_Sc)^0.5,
    May_k=(2841*Slope*MayV+2.02)*(600/May_Sc)^0.5,
    Jun_k=(2841*Slope*JunV+2.02)*(600/Jun_Sc)^0.5,
    Jul_k=(2841*Slope*JulV+2.02)*(600/Jul_Sc)^0.5,
    Aug_k=(2841*Slope*AugV+2.02)*(600/Aug_Sc)^0.5,
    Sep_k=(2841*Slope*SepV+2.02)*(600/Sep_Sc)^0.5,
    Oct_k=(2841*Slope*OctV+2.02)*(600/Oct_Sc)^0.5,
    Nov_k=(2841*Slope*NovV+2.02)*(600/Nov_Sc)^0.5,
    Dec_k=(2841*Slope*DecV+2.02)*(600/Dec_Sc)^0.5
  )

  df_2<-df[df$Slope>0.01,]
  df_2<-df_2%>%mutate(
    Jan_k=case_when(Jan_eD<=0.02~(2841*Slope*JanV+2.02)*(600/Jan_Sc)^0.5,
                    # Jan_eD<=0.02~(exp(0.35*log(Jan_eD)+3.10)*(600/Jan_Sc)^0.5),
                    Jan_eD>0.02~(exp(1.18*log(Jan_eD)+6.43)*(600/Jan_Sc)^0.5)),
    Feb_k=case_when(Feb_eD<=0.02~(2841*Slope*FebV+2.02)*(600/Feb_Sc)^0.5,
                    # Feb_eD<=0.02~(exp(0.35*log(Feb_eD)+3.10)*(600/Feb_Sc)^0.5),
                    Feb_eD>0.02~(exp(1.18*log(Feb_eD)+6.43)*(600/Feb_Sc)^0.5)),
    Mar_k=case_when(Mar_eD<=0.02~(2841*Slope*MarV+2.02)*(600/Mar_Sc)^0.5,
                    # Mar_eD<=0.02~(exp(0.35*log(Mar_eD)+3.10)*(600/Mar_Sc)^0.5),
                    Mar_eD>0.02~(exp(1.18*log(Mar_eD)+6.43)*(600/Mar_Sc)^0.5)),
    Apr_k=case_when(Apr_eD<=0.02~(2841*Slope*AprV+2.02)*(600/Apr_Sc)^0.5,
                    # Apr_eD<=0.02~(exp(0.35*log(Apr_eD)+3.10)*(600/Apr_Sc)^0.5),
                    Apr_eD>0.02~(exp(1.18*log(Apr_eD)+6.43)*(600/Apr_Sc)^0.5)),
    May_k=case_when(May_eD<=0.02~(2841*Slope*MayV+2.02)*(600/May_Sc)^0.5,
                    # May_eD<=0.02~(exp(0.35*log(May_eD)+3.10)*(600/May_Sc)^0.5),
                    May_eD>0.02~(exp(1.18*log(May_eD)+6.43)*(600/May_Sc)^0.5)),
    Jun_k=case_when(Jun_eD<=0.02~(2841*Slope*JunV+2.02)*(600/Jun_Sc)^0.5,
                    # Jun_eD<=0.02~(exp(0.35*log(Jun_eD)+3.10)*(600/Jun_Sc)^0.5),
                    Jun_eD>0.02~(exp(1.18*log(Jun_eD)+6.43)*(600/Jun_Sc)^0.5)),
    Jul_k=case_when(Jul_eD<=0.02~(2841*Slope*JulV+2.02)*(600/Jul_Sc)^0.5,
                    # Jul_eD<=0.02~(exp(0.35*log(Jul_eD)+3.10)*(600/Jul_Sc)^0.5),
                    Jul_eD>0.02~(exp(1.18*log(Jul_eD)+6.43)*(600/Jul_Sc)^0.5)),
    Aug_k=case_when(Aug_eD<=0.02~(2841*Slope*AugV+2.02)*(600/Aug_Sc)^0.5,
                    # Aug_eD<=0.02~(exp(0.35*log(Aug_eD)+3.10)*(600/Aug_Sc)^0.5),
                    Aug_eD>0.02~(exp(1.18*log(Aug_eD)+6.43)*(600/Aug_Sc)^0.5)),
    Sep_k=case_when(Sep_eD<=0.02~(2841*Slope*SepV+2.02)*(600/Sep_Sc)^0.5,
                    # Sep_eD<=0.02~(exp(0.35*log(Sep_eD)+3.10)*(600/Sep_Sc)^0.5),
                    Sep_eD>0.02~(exp(1.18*log(Sep_eD)+6.43)*(600/Sep_Sc)^0.5)),
    Oct_k=case_when(Oct_eD<=0.02~(2841*Slope*OctV+2.02)*(600/Oct_Sc)^0.5,
                    # Oct_eD<=0.02~(exp(0.35*log(Oct_eD)+3.10)*(600/Oct_Sc)^0.5),
                    Oct_eD>0.02~(exp(1.18*log(Oct_eD)+6.43)*(600/Oct_Sc)^0.5)),
    Nov_k=case_when(Nov_eD<=0.02~(2841*Slope*NovV+2.02)*(600/Nov_Sc)^0.5,
                    # Nov_eD<=0.02~(exp(0.35*log(Nov_eD)+3.10)*(600/Nov_Sc)^0.5),
                    Nov_eD>0.02~(exp(1.18*log(Nov_eD)+6.43)*(600/Nov_Sc)^0.5)),
    Dec_k=case_when(Dec_eD<=0.02~(2841*Slope*DecV+2.02)*(600/Dec_Sc)^0.5,
                    # Dec_eD<=0.02~(exp(0.35*log(Dec_eD)+3.10)*(600/Dec_Sc)^0.5),
                    Dec_eD>0.02~(exp(1.18*log(Dec_eD)+6.43)*(600/Dec_Sc)^0.5))
  )
  
  df<-rbind(df_1,df_2)
  rm(df_1,df_2)

  #saving the results
  write_csv(df[,c('COMID',paste0(month.abb,'_k'))],
            paste0(wd,'/output/table/MeritHydro/k_',str_pad(i,width=2,pad=0),'.csv'))
}
