library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(car)
library(caret)
library(gstat)
library(sp)
library(stringr)
library(corrplot)
library(nlme)
library(randomForest)
library(caTools)
library(scales)
library(MASS)
library(ggpubr)
library(viridis)

# set working directory
dir = 'D:/Research/globalEmission'

co2meas<-read_csv(paste(dir, 'dataset','directCO2', 'co2MeasRiverMons_noRes_2partsTog.csv', sep = '/'))
co2meas$lnCO2<-log(co2meas$CO2_uatm)

co2ms<-co2meas%>%group_by(siteNo)%>%summarise(
  Reference=Reference[1],
  SO=SO[1],
  Lat=mean(Lat),
  Long=mean(Long),
  CO2=mean(CO2_uatm),
  lnCO2=log(CO2))

#join SO, SO0,openwaterbodies; SO255, smallest stream sites
SO<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/SO_earthEnv.csv'))
SO[SO$SO%in%c(0),]$SO<--1 #-1,openwater 
SO[SO$SO%in%c(255),]$SO<-0 #smallest streams
#part2
SO2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/siteSO.csv'))
SO2<-SO2[,c('siteNo','strmOrder')]
names(SO2)[2]<-'SO'
SO<-rbind(SO,SO2)
co2meas<-left_join(co2meas,SO,by='siteNo')
rm(SO,SO2)
#rm open water
co2meas<-co2meas[(co2meas$SO!=-1),] #204 records removed

#join watershed area
basinarea=read_csv(paste0(dir,'/output/table/co2MeasRiverSites/basinArea.csv')) #km2
basinarea2=read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/basinArea.csv'))
basinarea<-rbind(basinarea,basinarea2)
co2meas<-left_join(co2meas,basinarea,by='siteNo')
co2ms<-left_join(co2ms,basinarea,by='siteNo')
rm(basinarea,basinarea2)

#join pop density, persons per km2
popdens<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/popdens.csv'))
popdens$popdens<-NA
popdens$popdens<-rowMeans(popdens[,2:5])
popdens2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/popdens.csv'))
popdens2$popdens<-NA
popdens2$popdens<-rowMeans(popdens2[,2:5])
popdens<-rbind(popdens,popdens2)
co2meas<-left_join(co2meas,popdens[,c(1,6)],by='siteNo')
co2ms<-left_join(co2ms,popdens[,c(1,6)],by='siteNo')
rm(popdens,popdens2)

#join precipitation and temperature
prectemp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/tempPrep.csv'))
prectemp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/tempPrep.csv'))
prectemp<-rbind(prectemp,prectemp2)
co2meas<-left_join(co2meas,prectemp,by='siteNo')
co2ms<-left_join(co2ms,prectemp,by='siteNo')
rm(prectemp,prectemp2)
co2meas<-
  co2meas%>%mutate(
  tavg_00=case_when(Mon%in%c('Jan')~tavg_01, #degree celcius
                    Mon%in%c('Feb')~tavg_02,
                    Mon%in%c('Mar')~tavg_03,
                    Mon%in%c('Apr')~tavg_04,
                    Mon%in%c('May')~tavg_05,
                    Mon%in%c('Jun')~tavg_06,
                    Mon%in%c('Jul')~tavg_07,
                    Mon%in%c('Aug')~tavg_08,
                    Mon%in%c('Sep')~tavg_09,
                    Mon%in%c('Oct')~tavg_10,
                    Mon%in%c('Nov')~tavg_11,
                    Mon%in%c('Dec')~tavg_12,
                    Mon%in%c('Ann')~tavg_ann),
  prec_00=case_when(Mon%in%c('Jan')~prec_01/31*365, #mm/mon
                    Mon%in%c('Feb')~prec_02/28*365,
                    Mon%in%c('Mar')~prec_03/31*365,
                    Mon%in%c('Apr')~prec_04/30*365,
                    Mon%in%c('May')~prec_05/31*365,
                    Mon%in%c('Jun')~prec_06/30*365,
                    Mon%in%c('Jul')~prec_07/31*365,
                    Mon%in%c('Aug')~prec_08/31*365,
                    Mon%in%c('Sep')~prec_09/30*365,
                    Mon%in%c('Oct')~prec_10/31*365,
                    Mon%in%c('Nov')~prec_11/30*365,
                    Mon%in%c('Dec')~prec_12/31*365,
                    Mon%in%c('Ann')~prec_ann) #mm/year
)
co2meas<-co2meas[,-which(names(co2meas)%in%paste0('tavg_',str_pad(c(1:12),2,pad=0)))]
co2meas<-co2meas[,-which(names(co2meas)%in%paste0('prec_',str_pad(c(1:12),2,pad=0)))]
co2ms<-co2ms[,-which(names(co2ms)%in%paste0('tavg_',str_pad(c(1:12),2,pad=0)))]
co2ms<-co2ms[,-which(names(co2ms)%in%paste0('prec_',str_pad(c(1:12),2,pad=0)))]

#join annual gpp, npp
anngppnpp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/annual_gpp_npp.csv'))
anngppnpp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/annual_gpp_npp.csv'))
anngppnpp<-rbind(anngppnpp,anngppnpp2)
co2meas<-left_join(co2meas,anngppnpp,by='siteNo')
co2ms<-left_join(co2ms,anngppnpp,by='siteNo')
rm(anngppnpp,anngppnpp2)

#join monthly_gpp_npp
mongppnpp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/monthly_gpp_npp.csv'))
mongppnpp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/monthly_gpp_npp.csv'))
mongppnpp<-rbind(mongppnpp,mongppnpp2)
co2meas<-left_join(co2meas,mongppnpp,by='siteNo')
rm(mongppnpp,mongppnpp2)
co2meas<-
  co2meas%>%mutate(
    GPP_00=case_when(Mon%in%c('Jan')~gpp_01,
                     Mon%in%c('Feb')~gpp_02,
                     Mon%in%c('Mar')~gpp_03,
                     Mon%in%c('Apr')~gpp_04,
                     Mon%in%c('May')~gpp_05,
                     Mon%in%c('Jun')~gpp_06,
                     Mon%in%c('Jul')~gpp_07,
                     Mon%in%c('Aug')~gpp_08,
                     Mon%in%c('Sep')~gpp_09,
                     Mon%in%c('Oct')~gpp_10,
                     Mon%in%c('Nov')~gpp_11,
                     Mon%in%c('Dec')~gpp_12,
                     Mon%in%c('Ann')~GPP),
    NPP_00=case_when(Mon%in%c('Jan')~npp_01,
                     Mon%in%c('Feb')~npp_02,
                     Mon%in%c('Mar')~npp_03,
                     Mon%in%c('Apr')~npp_04,
                     Mon%in%c('May')~npp_05,
                     Mon%in%c('Jun')~npp_06,
                     Mon%in%c('Jul')~npp_07,
                     Mon%in%c('Aug')~npp_08,
                     Mon%in%c('Sep')~npp_09,
                     Mon%in%c('Oct')~npp_10,
                     Mon%in%c('Nov')~npp_11,
                     Mon%in%c('Dec')~npp_12,
                     Mon%in%c('Ann')~NPP)
  )
co2meas<-co2meas[,-which(names(co2meas) %in% paste0('gpp_',str_pad(c(1:12),2,pad=0)))]
co2meas<-co2meas[,-which(names(co2meas) %in% paste0('npp_',str_pad(c(1:12),2,pad=0)))]

#join soil respiration
soilResp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/soilResp.csv'))
soilResp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/soilResp.csv'))
soilResp<-rbind(soilResp,soilResp2)
co2meas<-left_join(co2meas,soilResp,by='siteNo')
co2ms<-left_join(co2ms,soilResp,by='siteNo')
rm(soilResp,soilResp2)
co2meas<-
  co2meas%>%mutate(
    SR_00=case_when(Mon%in%c('Jan')~SR_01*365, #g C m-2 d-1
                      Mon%in%c('Feb')~SR_02*365,
                      Mon%in%c('Mar')~SR_03*365,
                      Mon%in%c('Apr')~SR_04*365,
                      Mon%in%c('May')~SR_05*365,
                      Mon%in%c('Jun')~SR_06*365,
                      Mon%in%c('Jul')~SR_07*365,
                      Mon%in%c('Aug')~SR_08*365,
                      Mon%in%c('Sep')~SR_09*365,
                      Mon%in%c('Oct')~SR_10*365,
                      Mon%in%c('Nov')~SR_11*365,
                      Mon%in%c('Dec')~SR_12*365,
                      Mon%in%c('Ann')~SR_ann) # g C m-2 yr-1
  )
co2meas<-co2meas[,-which(names(co2meas)%in%paste0('SR_',str_pad(c(1:12),2,pad=0)))]
co2ms<-co2ms[,-which(names(co2ms)%in%paste0('SR_',str_pad(c(1:12),2,pad=0)))]

#join wetland
wetland<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/wetland.csv'))
wetland2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/wetland.csv'))
wetland<-rbind(wetland,wetland2)
wetland$wetland<-wetland$wetland*100
co2meas<-left_join(co2meas,wetland,by='siteNo')
co2ms<-left_join(co2ms,wetland,by='siteNo')
rm(wetland,wetland2)

#join elevation (m) and slope (degree)
elevSlope<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/elevSlope.csv'))
elevSlope2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/elevSlope.csv'))
elevSlope<-rbind(elevSlope,elevSlope2)
elevSlope$slop<-tan(elevSlope$slop*pi/180) #m/m
co2meas<-left_join(co2meas,elevSlope,by='siteNo')
co2ms<-left_join(co2ms,elevSlope,by='siteNo')
rm(elevSlope,elevSlope2)

#join soil attributes
soilatts<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/soilatts.csv'))
soilatts2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/soilatts.csv'))
soilatts<-rbind(soilatts,soilatts2)
#GRAVEL %v, Sand %w, Silt %w, Clay %w, Density kg/L
#SOC %w, CEC cmol/kg, BS (Base Saturation) %, TEB (Total Exchangable Bases) cmol/kg
#CACO3 Calcium Carbonate, %w, CASO4 Gypsum, %w, ESP Sodicity or exchangable Na/CEC %
#ECE electrical conductivity, dS/m
soilatts<-soilatts[,c('siteNo','T_GRAVEL','S_SAND','T_SILT','T_CLAY','T_BULK_DENSITY',
                      'T_OC','T_PH_H2O', 'T_CEC_SOIL','T_BS','T_TEB',
                      'T_CACO3','T_CASO4','T_ESP','T_ECE')]
names(soilatts)<-c('siteNo','SGRAVEL','SSAND','SSILT','SCLAY','SDENSITY',
                   'SOC','SpH','SCEC','SBS','STEB','SCACO3','SCASO4',
                   'SSODICITY','SCONDUCTIVITY')
for(col in names(soilatts)){
  print(col);
  soilatts[soilatts[[col]]%in%c(0),col]<-NA
}
co2meas<-left_join(co2meas,soilatts,by='siteNo')
co2ms<-left_join(co2ms,soilatts,by='siteNo')
rm(soilatts,soilatts2)

ds<-co2meas
#do necessary transformation
ds$area<-log(ds$area)
#elev
ds$elev<-log(ds$elev)
#slop
ds$slop<-log(ds$slop)
#SOC
ds$SOC<-log(ds$SOC)
#SCEC
ds$SCEC<-log(ds$SCEC)
#STEB
ds$STEB<-log(ds$STEB)
#SCACO3
ds$SCACO3<-log(ds$SCACO3)
#SCASO4
ds$SCASO4<-log(ds$SCASO4)
#popdens
# ds[ds$popdens%in%c(0),]$popdens<-NA
ds$popdens<-log(ds$popdens)

#tranform for the annual dataset
co2ms$area<-log(co2ms$area)
co2ms$elev<-log(co2ms$elev)
co2ms$slop<-log(co2ms$slop)
#SOC, 8 zeros
co2ms$SOC<-log(co2ms$SOC)
#SCEC 8 zeros
co2ms$SCEC<-log(co2ms$SCEC)
#STEB,8 zeros
co2ms$STEB<-log(co2ms$STEB)
#SCACO3 543 zeros
co2ms$SCACO3<-log(co2ms$SCACO3)
#SCASO4 1000 zeros
co2ms$SCASO4<-log(co2ms$SCASO4)
#popdens 60 zeros
# co2ms[co2ms$popdens%in%c(0),]$popdens<-NA
co2ms$popdens<-log(co2ms$popdens)

####scatter plots####
cols_to_include<-c('lnCO2','tavg_ann','prec_ann','GPP','NPP','SR_ann',
                   'AR_ann','HR_ann','wetland','popdens','area','slop','elev','SGRAVEL',
                   'SSAND','SSILT','SCLAY','SDENSITY','SOC','SpH','SCEC',
                   'SBS','STEB','SCACO3','SCASO4','SSODICITY','SCONDUCTIVITY')
dsv<-co2ms[co2ms$popdens<log(300),cols_to_include]%>%
  gather(key='predictor',value=value,-lnCO2)
dsv<-dsv[!(is.na(dsv$value)),]
dsv<-dsv[!(is.infinite(dsv$value)),]
#rm extreme values
dsv<-dsv[!(dsv$predictor=='SpH'&dsv$value<3.7),] #rm SpH<3.7
dsv<-dsv[!(dsv$predictor=='SDENSITY'&dsv$value<1),] #rm SDENSITY<1
dsv<-dsv[!(dsv$predictor=='SSODICITY'&dsv$value>5),] #rm SSODICITY>5
dsv<-dsv[!(dsv$predictor=='SCONDUCTIVITY'&dsv$value>0.3),] #rm SCONDUCTIVITY>0.3

dsv$predictor<-
  factor(dsv$predictor,
         levels=c('tavg_ann','prec_ann',
                  'GPP','NPP','SR_ann','AR_ann','HR_ann',
                  'area','slop','elev','wetland','popdens',
                  'SGRAVEL','SSAND','SSILT','SCLAY','SDENSITY',
                  'SOC','SpH','SCEC','SBS','STEB','SCACO3','SCASO4',
                  'SSODICITY','SCONDUCTIVITY'),
         labels = c('Temperature (°C)','Precipitation (mm/yr)',
                    'GPP','NPP','Soil resp.','Autotrophic soil resp.','Heterotrophic soil resp.',
                    'Ln Watershed area (km2)','Ln Slope (unitless)','Ln Elevation (m)','Wetland (%)','Ln Pop dens. (people/km)',
                    'Soil gravel (%v)','Soil sand (%w)','Soil silt (%w)','Soil clay (%w)','Soil density (kg/m3)',
                    'Ln SOC (%w)','Soil pH','Ln Soil CEC (cmol/kg)','Base saturation (%)','Ln Soil TEB (cmol/kg)',
                    'Ln Calcium carbonate (%w)','Ln Soil gypsum (%w)',
                    'Soil sodiciy (%)','Soil conductivity (dS/m)'))

lm_r2<-
  function(predictor){
    fit=lm(lnCO2~value,dsv[(dsv$predictor==predictor)&(!is.na(dsv$predictor)),])
    fitsum<-summary(fit)
    r2=fitsum$adj.r.squared
    return(r2)
  }

pltText<-data.frame(R2=matrix(NA,26,1))
for(i in 1:26){
  print(i)
  r2=lm_r2(unique(dsv$predictor)[i])
  pltText[i,c('predictor')]<-unique(dsv$predictor)[i]
  pltText[i,c('R2')]<-round(r2,2)
  pltText[i,c('R2_Label')]<-
    as.character(as.expression(substitute(R^2*' = '*Rsqr, list(Rsqr=round(r2,2)))))
}

dsv%>%
  ggplot(aes(value,lnCO2))+
  geom_point(alpha=0.2,size=2)+
  geom_smooth(method=lm)+
  scale_y_continuous(limits=c(3,10),breaks=seq(4,10,by=2))+
  labs(y=expression('Ln '*CO[2]*' ('*mu*'atm)'))+
  geom_label(data=pltText,
            mapping=aes(x=Inf,y=3.8,label=R2_Label),
            hjust=1.2, parse=TRUE)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        panel.border=element_rect(colour = "black",fill=NA,size=1),
        strip.background=element_blank(),
        strip.text=element_text(size=12,color='red'))+
  facet_wrap(.~predictor,ncol=5,scales='free_x')
ggsave(paste0(dir,'/output/figure/CO2/linearRegressions_mon.jpg'),
       width=30,height=36,units='cm',device='jpeg')

####step-wise regression####
#1
ds[is.infinite(ds$popdens),]$popdens<-(-15)
m_min<-lm(lnCO2~1,data=ds)
m_max<-lm(lnCO2~tavg_00+prec_00+GPP_00+NPP_00+SR_00+
            slop+elev+wetland+popdens+
            SpH+SBS,data=ds)
m_red<-step(m_max,scope=c(upper=m_max,lower=m_min),
            trace=1)
summary(m_red) #0.4237
vif(m_red)
print(paste("AIC =",round(AIC(m_red),0)))
print(paste("N =",length(m_red$residuals)))

#2
m_min<-lm(lnCO2~1,data=ds)
m_max<-lm(lnCO2~tavg_00+prec_00+SR_00+
            slop+elev+wetland+popdens+
            SpH+SBS,data=ds)
m_red<-step(m_max,scope=c(upper=m_max,lower=m_min),
            trace=1)
summary(m_red) #0.4192
vif(m_red)
print(paste("AIC =",round(AIC(m_red),0)))

#3
m_min<-lm(lnCO2~1,data=ds)
m_max<-lm(lnCO2~tavg_00+prec_00+SR_00+
            slop+elev+wetland+popdens+
            SpH,data=ds)
m_red<-step(m_max,scope=c(upper=m_max,lower=m_min),
            trace=1)
summary(m_red) #0.4191
vif(m_red)
print(paste("AIC =",round(AIC(m_red),0)))

#4
m_min<-lm(lnCO2~1,data=ds)
m_max<-lm(lnCO2~prec_00+SR_00+
            slop+elev+wetland+popdens+
            SpH,data=ds)
m_red<-step(m_max,scope=c(upper=m_max,lower=m_min),
            trace=1)
summary(m_red) #0.4071
vif(m_red)
print(paste("AIC =",round(AIC(m_red),0)))

####Random Forest Model####
set.seed(0)
sa<-sample.split(1:5910,SplitRatio=0.75)
train<-subset(ds,sa)
test<-subset(ds,!sa)

#contruct the model
rfmod<-randomForest(lnCO2~tavg_00+prec_00+GPP_00+NPP_00+SR_00+area+
                      slop+elev+wetland+popdens+
                      SGRAVEL+SSAND+SOC+SpH+SBS,data=train,importance=TRUE,na.action=na.omit)

#model predictors importance
odr<-order(-rfmod$importance[,1])
importRF<-data.frame(rfmod$importance[odr,1])
importRF$pred<-row.names(importRF)
row.names(importRF)<-NULL
names(importRF)[1]<-'Import'
importRF$pred<-
  factor(importRF$pred,
         levels=c('slop','tavg_00','SR_00','elev','SpH','GPP_00',
                  'prec_00','NPP_00','SBS','SOC','popdens','SGRAVEL','area','SSAND','wetland'),
         labels=c('Slope','Temperature','Soil resp.','Elevation','Soil pH','GPP','Precipitation','NPP',
                  'Soil base saturation','SOC','Pop. dens.','Soil gravel',
                  'Watershed area','Soil sand','Wetland'))

importRF%>%ggplot(aes(x=pred,weight=Import))+
  geom_bar(fill='grey70')+
  labs(x='',y='Predictor importance')+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=1),
        axis.text.x=element_text(angle=90))

####model test and figure 1####
pred_rf<-predict(rfmod,newdata=test[,c('tavg_00','prec_00','GPP_00','NPP_00','SR_00','area','elev','slop',
                                       'wetland','popdens', 'SGRAVEL','SSAND','SOC','SpH','SBS')])
pred_lm<-predict(m_red,test[,c('SR_00','slop','elev','wetland','SpH','prec_00','popdens')])

mds<-data.frame(test$lnCO2,test$elev,test$Lat,pred_lm,pred_rf)
names(mds)<-c('lnco2','elev','Lat','pred_lm','pred_rf')
mds<-
  mds%>%mutate(
    res_lm=pred_lm-lnco2,
    res_rf=pred_rf-lnco2)
fit=lm(lnco2~pred_lm,mds)
fitsum=summary(fit)
fitsum#lm: 0.3814
fit=lm(lnco2~pred_rf,mds)
fitsum=summary(fit)
fitsum#rf: 0.755

fit=lm(exp(lnco2)~exp(pred_rf),mds)
fitsum=summary(fit)
fitsum#rf: 0.6767

pmain<-
  mds%>%
  ggplot(aes(x=exp(lnco2),y=exp(pred_rf)))+
  geom_point(size=4,alpha=1,aes(color=abs(Lat)))+
  geom_abline(intercept=0,slope=1,color='red',size=1)+
  annotate('text',500,11000,label=expression(R^2*' = 0.77'),size=6,color='red')+
  annotate('text',500,300,label=expression('1:1'),size=6,color='red')+
  labs(x=expression('Measured '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Predicted '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       color=expression('Abs. latitude (°)'))+
  scale_x_continuous(limits=c(200,20000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=2),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(200,20000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=2),
                     labels=trans_format('log10',math_format(10^.x)))+
  # scale_color_gradient(low='yellow2',high='green3')+
  scale_color_viridis(discrete=FALSE,direction=-1,alpha=0.5)+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        strip.background=element_blank(),
        legend.position=c(0.79,0.2),
        legend.key.size=unit(0.4,'cm'))

nfit<-fitdistr(mds$res_rf[!is.na(mds$res_rf)], "normal")
para<-nfit$estimate
perror<-
  mds%>%
  ggplot(aes(x=res_rf))+
  geom_histogram(bins=50,aes(y=..density..),color='grey50',size=0.5,fill='grey20')+
  stat_function(fun=dnorm,args=list(para[1],para[2]),color='red',size=1)+
  geom_segment(x=para[1],y=0,xend=para[1],yend=dnorm(para[1],para[1],para[2]),color='red',size=1)+
  geom_segment(x=para[1]-para[2],y=0,xend=para[1]-para[2],
               yend=dnorm(para[1]-para[2],para[1],para[2]),
               color='red',size=1,linetype=2)+
  geom_segment(x=para[1]+para[2],y=0,xend=para[1]+para[2],
               yend=dnorm(para[1]+para[2],para[1],para[2]),
               color='red',size=1,linetype=2)+
  annotate('text',-1.5,1.2,label=expression(mu*' = 0.009'),color='red',size=5)+
  annotate('text',-1.5,0.8,label=expression(delta*' = 0.42'),color='red',size=5)+
  annotate('text',para[1]+0.2,0.1+dnorm(para[1],para[1],para[2]),
           label=expression(mu),color='red',size=8)+
  annotate('text',para[1]+para[2]+0.4,dnorm(para[1]+para[2],para[1],para[2]),
           label=expression('+1'*delta),color='red',size=6)+
  annotate('text',para[1]-para[2]-0.4,dnorm(para[1]-para[2],para[1],para[2]),
           label=expression('-1'*delta),color='red',size=6)+
  labs(x=expression('Ln '*italic(p)*'C'*O[2]*' residuals ('*mu*'atm)'),
       y=expression('Density'))+
  scale_x_continuous(limits=c(-2,2))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        strip.background=element_blank())

co2meas<-co2meas%>%mutate(
  clmregion=case_when(abs(Lat)>=56~'Polar',
                      ((abs(Lat)<56)&abs(Lat)>=23.5)~'Temperate',
                      (abs(Lat)<23.5)~'Tropical'))

co2meas%>%group_by(clmregion)%>%summarise(
  co2=mean(CO2_uatm)
)

phist<-co2meas%>%
  ggplot(aes(CO2_uatm))+
  geom_histogram(binwidth=0.06,color='grey50',fill='grey20')+
  geom_segment(aes(x=390,y=0,xend=390,yend=150),linetype=2)+
  geom_text(aes(label=clmregion, x=30, y=310),
            data=data.frame(clmregion=c('Polar','Temperate','Tropical')),
            size=6,color='red',hjust='left') +
  labs(x=expression(italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Count'))+
  scale_x_continuous(limits=c(30,20000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(0,350))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text=element_blank())+
  facet_wrap(.~clmregion,ncol=1,strip.position='top')
ggsave(paste0(dir,'/output/figure/co2/fig1_mon.jpg'),
       ggarrange(phist,
                 ggarrange(pmain,perror,nrow=2,labels=c('B','C'),heights=c(2,1)),
                 ncol=2,widths=c(1,2),
                 labels='A'),
       width=20,height=16,units='cm',device='jpeg')

####comparing RF and RL errors on the testing dataset####
nfit_rf<-fitdistr(mds$res_rf[!is.na(mds$res_rf)], "normal")
para_rf<-nfit_rf$estimate

nfit_lm<-fitdistr(mds$res_lm[!is.na(mds$res_lm)], "normal")
para_lm<-nfit_lm$estimate

#p1
p1<-
  mds%>%
  ggplot(aes(x=exp(lnco2),y=exp(pred_lm)))+
  geom_point(size=4,alpha=0.6)+
  geom_smooth(method=lm,color='red',size=1)+
  annotate('text',1000,10000,label=expression(R^2*' = 0.38'),color='red',size=5)+
  annotate('text',10000,300,label=expression('RML model'),color='black',size=5)+
  labs(x=expression('Measured '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Predicted '*italic(p)*'C'*O[2]*' ('*mu*'atm)'))+
  scale_x_continuous(limits=c(200,25000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(200,25000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank())

p2<-mds%>%
  ggplot(aes(exp(lnco2),exp(pred_rf)))+
  geom_point(size=4,alpha=0.6)+
  geom_smooth(method=lm,color='red',size=1)+
  annotate('text',800,10000,label=expression(R^2*' = 0.77'),color='red',size=5)+
  annotate('text',10000,300,label=expression('RF model'),color='black',size=5)+
  labs(x=expression('Measured '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Predicted '*italic(p)*'C'*O[2]*' ('*mu*'atm)'))+
  scale_x_continuous(limits=c(200,25000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(200,25000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank())

p3<-
  mds%>%
  ggplot()+
  geom_histogram(aes(x=res_lm,y=..density..),color='grey50',bins=50)+
  stat_function(fun=dnorm,args=list(para_lm[1],para_lm[2]),color='red',size=1)+
  geom_segment(x=para_lm[1],y=0,xend=para_lm[1],yend=dnorm(para_lm[1],para_lm[1],para_lm[2]),color='red',size=1)+
  geom_segment(x=para_lm[1]-para_lm[2],y=0,xend=para_lm[1]-para_lm[2],
               yend=dnorm(para_lm[1]-para_lm[2],para_lm[1],para_lm[2]),
               color='red',size=1,linetype=2)+
  geom_segment(x=para_lm[1]+para_lm[2],y=0,xend=para_lm[1]+para_lm[2],
               yend=dnorm(para_lm[1]+para_lm[2],para_lm[1],para_lm[2]),
               color='red',size=1,linetype=2)+
  annotate('text',2,2,label=expression('RML model'),color='black',size=5)+
  annotate('text',-2,1.5,label=expression(mu*' = -0.02'),color='red',size=5)+
  annotate('text',-2,1.2,label=expression(delta*' = 0.69'),color='red',size=5)+
  annotate('text',para_lm[1]+0.4,0.1+dnorm(para_lm[1],para_lm[1],para_lm[2]),
           label=expression(mu),color='red',size=8)+
  annotate('text',para_lm[1]+para_lm[2]+0.8,dnorm(para_lm[1]+para_lm[2],para_lm[1],para_lm[2]),
           label=expression('+1'*delta),color='red',size=6)+
  annotate('text',para_lm[1]-para_lm[2]-0.8,dnorm(para_lm[1]-para_lm[2],para_lm[1],para_lm[2]),
           label=expression('-1'*delta),color='red',size=6)+
  scale_y_continuous(limits=c(0,2.3))+
  scale_x_continuous(limits=c(-3,3))+
  labs(x=expression('Ln '*italic(p)*'C'*O[2]*' residual'),
       y=expression('Density'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank())

p4<-
  mds%>%
  ggplot()+
  geom_histogram(aes(x=res_rf,y=..density..),color='grey50',bins=50)+
  stat_function(fun=dnorm,args=list(para_rf[1],para_rf[2]),color='red',size=1)+
  geom_segment(x=para_rf[1],y=0,xend=para_rf[1],yend=dnorm(para_rf[1],para_rf[1],para_rf[2]),color='red',size=1)+
  geom_segment(x=para_rf[1]-para_rf[2],y=0,xend=para_rf[1]-para_rf[2],
               yend=dnorm(para_rf[1]-para_rf[2],para_rf[1],para_rf[2]),
               color='red',size=1,linetype=2)+
  geom_segment(x=para_rf[1]+para_rf[2],y=0,xend=para_rf[1]+para_rf[2],
               yend=dnorm(para_rf[1]+para_rf[2],para_rf[1],para_rf[2]),
               color='red',size=1,linetype=2)+
  annotate('text',2,2,label=expression('RF model'),color='black',size=5)+
  annotate('text',-2,1.5,label=expression(mu*' = 0.009'),color='red',size=5)+
  annotate('text',-2,1.2,label=expression(delta*' = 0.42'),color='red',size=5)+
  annotate('text',para_rf[1]+0.4,0.1+dnorm(para_rf[1],para_rf[1],para_rf[2]),
           label=expression(mu),color='red',size=8)+
  annotate('text',para_rf[1]+para_rf[2]+0.8,dnorm(para_rf[1]+para_rf[2],para_rf[1],para_rf[2]),
           label=expression('+1'*delta),color='red',size=6)+
  annotate('text',para_rf[1]-para_rf[2]-0.8,dnorm(para_rf[1]-para_rf[2],para_rf[1],para_rf[2]),
           label=expression('-1'*delta),color='red',size=6)+
  scale_y_continuous(limits=c(0,2.3))+
  scale_x_continuous(limits=c(-3,3))+
  labs(x=expression('Ln '*italic(p)*'C'*O[2]*' residual'),
       y=expression('Density'))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank())

p5<-mds%>%
  ggplot(aes(exp(lnco2),res_lm))+
  geom_point(size=4,alpha=0.6)+
  geom_hline(yintercept=0,linetype=2,color='red',size=1)+
  annotate('text',300,-0.8,label='Zero line',size=5,color='red')+
  annotate('text',10000,2.5,label=expression('RML model'),color='black',size=5)+
  labs(x=expression(italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Model residual'))+
  scale_x_continuous(limits=c(100,30000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(-3,3),breaks=seq(-3,3,2))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank())

p6<-
  mds%>%
  ggplot(aes(x=exp(lnco2),y=res_rf))+
  geom_point(size=4,alpha=0.6)+
  geom_hline(yintercept=0,linetype=2,color='red',size=1)+
  annotate('text',300,-0.8,label='Zero line',size=5,color='red')+
  annotate('text',10000,2.5,label=expression('RF model'),color='black',size=5)+
  labs(x=expression(italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Model residual'))+
  scale_x_continuous(limits=c(100,30000),
                     trans=log_trans(base=10),
                     breaks=trans_breaks('log10',function(x) 10^x,n=3),
                     labels=trans_format('log10',math_format(10^.x)))+
  scale_y_continuous(limits=c(-3,3),breaks=seq(-3,3,2))+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        legend.position=c(0.2,0.8),
        legend.title=element_blank())
ggsave(paste0(dir,'/output/figure/co2/rF_lm_modelcomp_mon.jpg'),
       plot=grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2),
       width=20,height=25,units='cm',device='jpeg')

####10-fold cross validation####
set.seed(0)

flds<-createFolds(1:5910,k=10,list=TRUE,returnTrain = FALSE)
names(flds)

rst=data.frame(matrix(NA,nrow=10,ncol=2))
names(rst)<-c('r2','error')

for(i in 1:10){
  print(i)
  train<-ds[!(1:5910 %in% flds[[i]]),]
  test<-ds[flds[[i]],]
  
  #contruct the model
  rfmod<-randomForest(lnCO2~tavg_00+prec_00+GPP_00+NPP_00+SR_00+area+
                        slop+elev+wetland+popdens+
                        SGRAVEL+SSAND+SOC+SpH+SBS,data=train,importance=TRUE,na.action=na.omit,ntree=200)
  pred_rf<-predict(rfmod,newdata=test[,c('tavg_00','prec_00','GPP_00','NPP_00','SR_00','area','elev','slop',
                                         'wetland','popdens', 'SGRAVEL','SSAND','SOC','SpH','SBS')])
 
  mds<-data.frame(test$lnCO2,test$elev,test$Lat,pred_rf)
  names(mds)<-c('lnco2','elev','Lat','pred_rf')
  mds<-
    mds%>%mutate(
      res_rf=pred_rf-lnco2)
  fit=lm(lnco2~pred_rf,mds)
  fitsum=summary(fit)
  # print(paste0('R2=',fitsum$r.squared))
  rst[i,]$r2<-fitsum$r.squared
  nfit<-fitdistr(mds$res_rf[!is.na(mds$res_rf)], "normal")
  para<-nfit$estimate
  # print(para[2])
  rst[i,]$error<-para[2]
}

mean(rst$r2)
sd(rst$r2)
mean(rst$error)
sd(rst$error)