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
  
  #calc k ray12
  # print(sum(df$Slope>0.01)/nrow(df))
  # if(i==1){highSlopeFL<-sum(df$Slope>0.01)}else{highSlopeFL<-highSlopeFL+sum(df$Slope>0.01)}
  # if(i==1){totFL<-nrow(df)}else{totFL<-totFL+nrow(df)}
  
  # df<-df%>%mutate(
  #   Jan_k=(2841*Slope*JanV+2.02)*(600/Jan_Sc)^0.5,
  #   Feb_k=(2841*Slope*FebV+2.02)*(600/Feb_Sc)^0.5,
  #   Mar_k=(2841*Slope*MarV+2.02)*(600/Mar_Sc)^0.5,
  #   Apr_k=(2841*Slope*AprV+2.02)*(600/Apr_Sc)^0.5,
  #   May_k=(2841*Slope*MayV+2.02)*(600/May_Sc)^0.5,
  #   Jun_k=(2841*Slope*JunV+2.02)*(600/Jun_Sc)^0.5,
  #   Jul_k=(2841*Slope*JulV+2.02)*(600/Jul_Sc)^0.5,
  #   Aug_k=(2841*Slope*AugV+2.02)*(600/Aug_Sc)^0.5,
  #   Sep_k=(2841*Slope*SepV+2.02)*(600/Sep_Sc)^0.5,
  #   Oct_k=(2841*Slope*OctV+2.02)*(600/Oct_Sc)^0.5,
  #   Nov_k=(2841*Slope*NovV+2.02)*(600/Nov_Sc)^0.5,
  #   Dec_k=(2841*Slope*DecV+2.02)*(600/Dec_Sc)^0.5
  # )

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

  # df_2a<-df_2%>%mutate(
  #   Jan_k=(2841*Slope*JanV+2.02)*(600/Jan_Sc)^0.5,
  #   Feb_k=(2841*Slope*FebV+2.02)*(600/Feb_Sc)^0.5,
  #   Mar_k=(2841*Slope*MarV+2.02)*(600/Mar_Sc)^0.5,
  #   Apr_k=(2841*Slope*AprV+2.02)*(600/Apr_Sc)^0.5,
  #   May_k=(2841*Slope*MayV+2.02)*(600/May_Sc)^0.5,
  #   Jun_k=(2841*Slope*JunV+2.02)*(600/Jun_Sc)^0.5,
  #   Jul_k=(2841*Slope*JulV+2.02)*(600/Jul_Sc)^0.5,
  #   Aug_k=(2841*Slope*AugV+2.02)*(600/Aug_Sc)^0.5,
  #   Sep_k=(2841*Slope*SepV+2.02)*(600/Sep_Sc)^0.5,
  #   Oct_k=(2841*Slope*OctV+2.02)*(600/Oct_Sc)^0.5,
  #   Nov_k=(2841*Slope*NovV+2.02)*(600/Nov_Sc)^0.5,
  #   Dec_k=(2841*Slope*DecV+2.02)*(600/Dec_Sc)^0.5
  # )
  # 
  # if(i!=1){df2<-rbind(df2,df_2)}else{df2<-df_2}
  # if(i!=1){df2a<-rbind(df2a,df_2a)}else{df2a<-df_2a}
  
  df<-rbind(df_1,df_2)
  rm(df_1,df_2)

  #saving the results
  write_csv(df[,c('COMID',paste0(month.abb,'_k'))],
            paste0(wd,'/output/table/MeritHydro/k_',str_pad(i,width=2,pad=0),'.csv'))
}

print(highSlopeFL/totFL*100)
# #k increase on high-slope streams 30%, 18 m d-1 to 24 m d-1 23.68/18*100; 38%: 31.8/23.1 m d-1
# mean(df2[,paste0(month.abb,'_k')]%>%rowMeans(),na.rm=T)/mean(df2a[,paste0(month.abb,'_k')]%>%rowMeans(),na.rm=T)*100
# write_csv(df_2all[,c('COMID','Slope')],paste0(wd,'/output/figure/k/COMID_eD002plus.csv'))

####k spatial distribution####
#putting up the k dataset
for(i in 1:8){
  print(paste0('Processing ',i))
  df_i=read.dbf(paste0(wd,'/dataset/MeritHydroLinPan/pfaf_0',i,'_riv_3sMERIT.dbf'))
  df_i<-df_i[,c('COMID','strmOrder','Slope','Length')]
  k_i<-read_csv(paste0(wd,'/output/table/MeritHydro/k_',str_pad(i,width=2,pad=0),'.csv'))
  df_i<-left_join(df_i,k_i,by='COMID')
  rm(k_i)
  if(i==1){df<-df_i}else{df<-rbind(df,df_i)}
}
# write_csv(df[,c('COMID',paste0(month.abb,'_k'))],paste0(wd,'/output/table/MeritHydro/k.csv'))
df$yea_k<-df[,paste0(month.abb,'_k')]%>%rowMeans()

#North American k SO changes
df[df$COMID>10000000&df$COMID<=80000000,]%>%
  group_by(strmOrder)%>%summarise(
    yea_k=mean(yea_k,na.rm=TRUE)
  )%>%ggplot(aes(strmOrder,yea_k))+
  geom_point(size=4)+
  scale_x_continuous(breaks=seq(1,9))+
  scale_y_continuous(limits=c(0,11),breaks=seq(0,10,2))+
  labs(x=expression('Stream order'),
       y=expression('Average k (m '*d^-1*')'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA, size=0.6))
ggsave(paste0(wd,'/output/figure/k/SO_k.png'),
       width=12,height=9,units='cm')

# North America
k_na<-
read_csv('D:/research/completed writings/Journal Papers/Liu and Raymond 2018 LOL/hydro_response/output/results/lm_allQ.csv')

k_na<-k_na[!(k_na$SO==9&k_na$k>20),]
k_na<-k_na[!(k_na$SO==8&k_na$k>20),]
k_na<-k_na[!(k_na$SO==7&k_na$k>20),]

k_na_<-k_na%>%group_by(SO)%>%summarise(k=mean(k))
k_na_$SO<-k_na_$SO-2
k_na_<-k_na_[k_na_$SO>0,]

k_na_GRADES<-df[df$COMID>70000000&df$COMID<=80000000,]%>%
  group_by(strmOrder)%>%summarise(yea_k=mean(yea_k,na.rm=TRUE))

left_join(k_na_,k_na_GRADES,by=c('SO'='strmOrder'))%>%
  gather(key=tp,value=k,-SO)%>%
  ggplot(aes(SO,k,color=tp))+
  geom_point()
  # ggplot(aes(k,yea_k))+
  # geom_point()+
  # geom_abline(slope=1,intercept=1)

#literature versus modeled k
k_lit<-read_csv(paste0(wd,'/dataset/directCO2/co2MeasRiverMons_2partsTog.csv'))
k_lit<-k_lit[!is.na(k_lit$k600_md),]
k_lit<-k_lit[!k_lit$siteType%in%c('lake','reservoir'),]
#join stream order
SO_part1<-read_csv(paste0(wd,'/output/table/co2MeasRiverSites/SO_earthEnv.csv'))
comid_part1<-read_excel(paste0(wd,'/dataset/directCO2/sitesCOMIDmatch.xlsx'))
comid_part1<-comid_part1[,c('siteNo','COMID')]
SO_part1<-left_join(SO_part1,comid_part1,by='siteNo')
SO_part1<-SO_part1[,c('siteNo','COMID','SO')]
SO_part2<-read_csv(paste0(wd,'/output/table/co2MeasRiverSites/part2/siteSO.csv'))
names(SO_part2)<-c('siteNo','COMID','SO')
SO<-rbind(SO_part1,SO_part2)
rm(SO_part1,comid_part1,SO_part2)
k_lit<-left_join(k_lit,SO,by='siteNo')
k_lit<-k_lit[k_lit$SO<10,]

k_lit$COMID<-as.integer(k_lit$COMID)
#join model k
k_lit<-left_join(k_lit,df[,c('COMID','yea_k')],by='COMID')

table(k_lit$Reference)
mean(k_lit[k_lit$Reference%in%c('Borges2015NatureGeosci'),]$yea_k,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Alin2011JGR'),]$yea_k,na.rm=TRUE)
mean(k_lit[k_lit$Reference%in%c('Alin2011JGR')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Ran2017JGR'),]$yea_k,na.rm=TRUE)
mean(k_lit[k_lit$Reference%in%c('Ran2017JGR')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Liu2016GBC'),]$yea_k,na.rm=TRUE)
mean(k_lit[k_lit$Reference%in%c('Liu2016GBC')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Amaral2019Biogeochem')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Rocher-Ros2019LOL')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)

mean(k_lit[k_lit$Reference%in%c('Teodoru2015Biogeosci')&k_lit$k600_md<100,]$k600_md,na.rm=TRUE)


k_lit%>%group_by(Reference)%>%summarise(
  klit=mean(k600_md,na.rm=TRUE),
  kmod=mean(yea_k,na.rm=TRUE)
)%>%View()


k_lit<-k_lit[!k_lit$Reference%in%c('Rocher-Ros2019LOL'),]

k_lit[k_lit$Reference%in%c('Alin2011JGR'),]%>%group_by(SO)%>%summarise(
  k=mean(k600_md)
)%>%ggplot(aes(SO,k))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,10,1))

  ggplot(aes(SO,k600_md))+
  geom_point()
  
k_lit[k_lit$Reference%in%c('Teodoru2015Biogeosci'),]%>%group_by(SO)%>%summarise(
    k=mean(k600_md)
  )%>%ggplot(aes(SO,k))+
    geom_point()+
    scale_x_continuous(breaks = seq(0,10,1))

k_lit[k_lit$Reference%in%c('Ran2017JGR'),]%>%group_by(SO)%>%summarise(
  k=mean(k600_md)
)%>%ggplot(aes(SO,k))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,10,1))

k_lit[k_lit$Reference%in%c('Amaral2019Biogeochem'),]%>%group_by(SO)%>%summarise(
  k=mean(k600_md)
)%>%ggplot(aes(SO,k))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,10,1))



k_lit%>%group_by(SO)%>%summarise(
  k=mean(k600_md)
)%>%ggplot(aes(SO,k))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,10,1))

table(k_lit$Reference)

k_lit%>%group_by(SO)%>%summarise(
  k=mean(k600_md)
)%>%ggplot(aes(SO,k))+
  geom_point()+
  scale_x_continuous(breaks = seq(0,10,1))

mean(df[df$strmOrder>=4,]$yea_k) #3.3 m d-1
mean(df[df$strmOrder<4,]$yea_k,na.rm=TRUE) #8.2



df%>%group_by(strmOrder)%>%summarise(
  yea_k=mean(yea_k,na.rm=TRUE)
)%>%ggplot(aes(strmOrder,yea_k))+
  geom_point(size=4)+
  scale_x_continuous(breaks=seq(1,9))+
  labs(x=expression('Stream order'),
       y=expression('Average k (m '*d^-1*')'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA, size=0.6))
ggsave(paste0(wd,'/output/figure/k/SO_k.png'),
       width=12,height=9,units='cm')

df[df$strmOrder>=3,]$yea_k%>%mean(na.rm=TRUE)

####comparing NHDplus Slope and GRADES river slope at UGSG gages####
#nhd stream order
nhd<-read_csv(paste0(wd,'/dataset/NHDplusdata/NHDFlowline_Network_1.csv'))
nhd<-nhd[,c('COMID','StreamOrde','FCODE')]
nhd<-nhd[nhd$FCODE%in%c(46000,46003,46006),]
#nhd slope
nhd5<-read_csv(paste0(wd,'/dataset/NHDplusdata/NHDFlowline_Network_5.csv'))
nhd5<-nhd5[,c('COMID','SLOPE')]
nhd5<-nhd5[nhd5$SLOPE>0,]
nhd<-inner_join(nhd,nhd5,by='COMID')
rm(nhd5)
names(nhd)<-c('nhdCOMID','StreamOrde','FCODE','SLOPE')
#comids
gauges<-read_csv(paste0(wd,'/output/table/usgsGaugeMeritFL/usgsStreamRiverSites_sjoined.csv'))
gauges<-gauges[,c('FLComID','COMID')]
names(gauges)<-c('nhdCOMID','meritCOMID')
#adding meritCOMID
slopCOMP<-inner_join(gauges,nhd,by='nhdCOMID')
#join merit Slope
df=read.dbf(paste0(wd,'/dataset/MeritHydroLinPan/pfaf_0',7,'_riv_3sMERIT.dbf'))
slopCOMP<-inner_join(slopCOMP,df[,c('COMID','Slope')],by=c('meritCOMID'='COMID'))
slopCOMP<-slopCOMP[slopCOMP$SLOPE!=0.00001,]
slopCOMP<-slopCOMP[slopCOMP$Slope>=0.00001,]
#tidy the dataset
slopCOMP_v<-slopCOMP%>%gather(key='Method',value='Slope',-nhdCOMID,-meritCOMID,-StreamOrde,-FCODE)
slopCOMP_v$Method<-
  factor(slopCOMP_v$Method,labels=c('GRADES','NHDplus'))

median(slopCOMP$SLOPE)
median(slopCOMP$SLOPE)

p1<-slopCOMP%>%
  ggplot(aes(SLOPE,Slope))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0,size=1,color='red')+
  scale_x_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x),
                     labels=trans_format('log10',math_format(10^.x)))+
  annotate('text',x=0.03,y=0.00003,label='1:1 Line',color='red',size=5)+
  labs(x='NHDplus slope',y='GRADES slope')+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,color='black',size=0.5))
p2<-slopCOMP_v%>%
  ggplot(aes(Slope,fill=Method))+
  geom_density(alpha=0.4)+
  scale_x_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x),
                     labels=trans_format('log10',math_format(10^.x)))+
  labs(x=expression('Slope'),fill=expression(''),
       y=expression('Density'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,color='black',size=0.5),
        legend.position=c(0.2,0.8))
ggsave(paste0(wd,'/output/figure/k/comparing slope.png'), grid.arrange(p1,p2,nrow=1),  width=24,height=9,units='cm')

df_2[df_2$Jan_eD<=0.02,]%>%ggplot(aes(Jan_k,Jan_k_U))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1,intercept=0)+
  labs(x=expression('k predicted from Raymond et al. (m '*d^-1*')'),
       y=expression('k predicted from Ulseth et al. (m '*d^-1*')'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.6),
        legend.position=c(0.2,0.8))
ggsave(paste0(wd,'/output/figure/k/comparing_k_b0.02.png'),
       width=12,height=9,units='cm')

#Ulseth19
ulseth<-read_excel(paste0(wd,'/dataset/Ulseth19_gasTransferMontain/41561_2019_324_MOESM2_ESM.xlsx'),sheet=2)

quantile(ulseth[ulseth$eD>0.02,]$slope,0.123)
quantile(ulseth[ulseth$eD<=0.02,]$slope,0.88)

ulseth%>%
  ggplot(aes(slope,fill=eD>0.02))+
  geom_density(color=NA,alpha=0.5)+
  geom_vline(xintercept=0.01,linetype=2,color='red')+
  annotate('text',x=0.035,y=0.7,label='Slope: 0.01',color='red',size=4)+
  scale_x_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x),
                     labels=trans_format('log10',math_format(10^.x)))+
  labs(x=expression('Slope'),
       y=expression('Density'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.6),
        legend.position=c(0.2,0.8))
ggsave(paste0(wd,'/output/figure/k/ulsethSlopes.png'),
       width=12,height=9,units='cm')

####k is sensative to slope cutoff####
slopSensT<-read_excel(paste0(wd,'/output/table/regionSurfArea/basinSurfaceAreaSummary.xlsx'),
                      sheet='slopeCutoffSensTest')
names(slopSensT)<-c('cutoff_1','cutoff_2','k','percReaches')
slopSensT<-slopSensT[slopSensT$cutoff_1>0,]

p1<-slopSensT%>%ggplot(aes(cutoff_1,percReaches))+
  geom_point(size=4)+
  labs(x='Slope Cutoff',
       y=expression('% GRADES Reaches as High Slope Streams'))+
  scale_x_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x, n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.6),
        legend.position=c(0.2,0.8))

p2<-slopSensT%>%ggplot(aes(cutoff_1,k))+
  geom_point(size=4)+
  labs(x='Slope Cutoff',
       y=expression('k (m '*d^-1*')'))+
  scale_x_continuous(trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x, n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(8,12),
                     breaks=seq(8,12,2))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.6),
        legend.position=c(0.2,0.8))
p<-grid.arrange(p1,p2,ncol=2)
ggsave(paste0(wd,'/output/figure/k/slopSensT.png'),
       plot=p,
       width=24,height=9,units='cm')

