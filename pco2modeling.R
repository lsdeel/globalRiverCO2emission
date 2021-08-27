library(readr)
library(readxl)
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
rm(popdens,popdens2)

#join precipitation and temperature
prectemp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/tempPrep.csv'))
prectemp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/tempPrep.csv'))
prectemp<-rbind(prectemp,prectemp2)
co2meas<-left_join(co2meas,prectemp,by='siteNo')
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

#join annual gpp, npp
anngppnpp<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/annual_gpp_npp.csv'))
anngppnpp2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/annual_gpp_npp.csv'))
anngppnpp<-rbind(anngppnpp,anngppnpp2)
co2meas<-left_join(co2meas,anngppnpp,by='siteNo')
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

#join wetland
wetland<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/wetland.csv'))
wetland2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/wetland.csv'))
wetland<-rbind(wetland,wetland2)
wetland$wetland<-wetland$wetland*100
co2meas<-left_join(co2meas,wetland,by='siteNo')
rm(wetland,wetland2)

#join elevation (m) and slope (degree)
elevSlope<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/elevSlope.csv'))
elevSlope2<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/part2/elevSlope.csv'))
elevSlope<-rbind(elevSlope,elevSlope2)
elevSlope$slop<-tan(elevSlope$slop*pi/180) #m/m
co2meas<-left_join(co2meas,elevSlope,by='siteNo')
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
rm(soilatts,soilatts2)

#join flow slope and Q
flslopeQ<-read_csv(paste0(dir,'/output/table/co2MeasRiverSites/FLslopeQ.csv'))
flslopeQ$V<-exp(-1.06+0.12*log(flslopeQ$annQ))
flslopeQ$k<-2841*flslopeQ$Slope*flslopeQ$V+2.02
co2meas<-left_join(co2meas,flslopeQ[c('siteNo','k')],by='siteNo')
rm(flslopeQ)

# remove influencial sites
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

####histograms for CO2 and predictors####
cols_to_include<-c('lnCO2','area','popdens','tavg_00','prec_00','GPP_00','NPP_00',
                   'AR_ann','HR_ann','SR_00','wetland','elev','slop','SGRAVEL',
                   'SSAND','SSILT','SCLAY','SDENSITY','SOC','SpH','SCEC','SBS',
                   'STEB','SCACO3','SCASO4','SSODICITY','SCONDUCTIVITY')
dsv<-ds[,which(names(ds) %in% cols_to_include)]
dsv<-dsv%>%gather(key='predictor',value=value)
dsv<-dsv[!is.na(dsv$value),]

#rm extreme values
dsv<-dsv[!(dsv$predictor=='SpH'&dsv$value<quantile(ds$SpH,0.01,na.rm=TRUE)),] #rm SpH<3.7
dsv<-dsv[!(dsv$predictor=='SDENSITY'&dsv$value<1),] #rm SDENSITY<1
dsv<-dsv[!(dsv$predictor=='SSODICITY'&dsv$value>5),] #rm SSODICITY>5
dsv<-dsv[!(dsv$predictor=='SCONDUCTIVITY'&dsv$value>0.3),] #rm SCONDUCTIVITY>0.3

#format predictors
dsv$predictor<-
  factor(dsv$predictor,
         levels=c('lnCO2','tavg_00','prec_00',
                  'GPP_00','NPP_00','SR_00','AR_ann','HR_ann',
                  'area','slop','elev','wetland','popdens',
                  'SGRAVEL','SSAND','SSILT','SCLAY','SDENSITY',
                  'SOC','SpH','SCEC','SBS','STEB','SCACO3','SCASO4',
                  'SSODICITY','SCONDUCTIVITY'),
         labels = c('Ln CO2','Temperature (°C)','Precipitation (mm/yr)',
                    'GPP','NPP','Soil resp.','Autotrophic soil resp.','Heterotrophic soil resp.',
                    'Ln Watershed area (km2)','Ln Slope (unitless)','Ln Elevation (m)','Wetland (%)','Ln pop. dens. (people/km2)',
                    'Soil gravel (%v)','Soil sand (%w)','Soil silt (%w)','Soil clay (%w)','Soil density (kg/m3)',
                    'Ln SOC (%w)','Soil pH','Ln Soil CEC (cmol/kg)','Base saturation (%)','Ln Exch. bases (cmol/kg)',
                    'Ln Calcium carbonate (%w)','Ln Soil gypsum (%w)',
                    'Soil sodicity (%)','Soil conductivity (dS/m)'))

#histograms
dsv%>%ggplot(aes(value))+
  geom_histogram(bins=50)+
  theme_classic()+
  theme(
    axis.title.x=element_blank(),
    panel.border=element_rect(colour="black",fill=NA,size=1),
    strip.background=element_blank(),
    strip.text=element_text(size=12,color='red')
  )+
  facet_wrap(.~predictor,scales='free')
# ggsave(paste0(dir,'/output/figure/co2/hist_pred_mon.png'),
#        width=30, height=20, units='cm')


####scatter plots####
cols_to_include<-c('lnCO2','tavg_00','prec_00','GPP_00','NPP_00','SR_00',
                   'AR_ann','HR_ann','wetland','popdens','area','slop','elev','SGRAVEL',
                   'SSAND','SSILT','SCLAY','SDENSITY','SOC','SpH','SCEC',
                   'SBS','STEB','SCACO3','SCASO4','SSODICITY','SCONDUCTIVITY')

dsv<-ds[ds$popdens<log(300),cols_to_include]%>%
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
         levels=c('tavg_00','prec_00',
                  'GPP_00','NPP_00','SR_00','AR_ann','HR_ann',
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
            mapping=aes(x=Inf,y=3.5,label=R2_Label),
            hjust=1.2, parse=TRUE)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        panel.border=element_rect(colour = "black",fill=NA,size=1),
        strip.background=element_blank(),
        strip.text=element_text(size=12,color='red'))+
  facet_wrap(.~predictor,ncol=5,scales='free_x')
# ggsave(paste0(dir,'/output/figure/co2/linearRegressions_monMean.jpg'),
#        width=30,height=36,units='cm',device='jpeg')


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

####lm model validation####
realpred<-
data.frame(modelvalue=m_red[["fitted.values"]],
           modelres=m_red[["residuals"]],
           realvalue=m_red[["model"]][["lnCO2"]],
           soilResp=m_red[["model"]][["SR_00"]],
           slop=m_red[["model"]][["slop"]],
           elev=m_red[["model"]][["elev"]],
           wetland=m_red[["model"]][["wetland"]],
           SpH=m_red[["model"]][["SpH"]],
           ds[!(1:5910%in%unname(unclass(m_red$na.action))), 
              c('siteNo','SO','Lat','area','popdens','tavg_ann','prec_ann','GPP','NPP','AR_ann','HR_ann','SGRAVEL', 'SSAND','SBS')])

summary(lm(modelvalue~realvalue,data=realpred))#0.4075

#predictors versus residuals
realpredv<-
  realpred[,c('modelres','soilResp','slop','elev','wetland','SpH')]%>%
  gather(key=predictor,value=value,-modelres)

realpredv$predictor<-
  factor(realpredv$predictor,
         levels=c('soilResp','slop','elev','wetland','SpH'),
         labels=c('Soil resp.','Ln Slope (unitless)',
                  'Ln Elevation (m)', 'Wetland (%)','Soil pH'))
realpredv<-realpredv[!(realpredv$predictor=='Soil pH' & realpredv$value<4),]

realpredv%>%
  ggplot(aes(value,modelres))+
  geom_point(alpha=0.4,size=1)+
  geom_hline(yintercept=0,color='blue')+
  labs(y=expression('Residuals'))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  facet_wrap(.~predictor,ncol=3,scales='free_x')
# ggsave(paste0(dir,'/output/figure/co2/residual_predictors_mon.png'),
#        width=24,height=12,units='cm',device='png')

#residuals versus non-predictors
realpredv<-
  realpred[,c('modelres','Lat','area','popdens','tavg_ann','prec_ann','GPP','NPP',
              'AR_ann','HR_ann','SGRAVEL','SSAND','SBS')]
# realpredv$popdens<-log(realpredv$popdens)
realpredv<-realpredv%>%gather(key=predictor,value=value,-modelres)

realpredv$predictor<-
  factor(realpredv$predictor,
         levels=c('Lat','area','popdens','tavg_ann','prec_ann','GPP','NPP',
                  'AR_ann','HR_ann','SGRAVEL','SSAND','SBS'),
         labels=c('Latitude','Watershed area','Pop density','Temperature','Precipitation','GPP','NPP',
                  'Autotrophic resp.', 'Heterotrophic resp.',
                  'Soil gravel (%v)','Soil sand (%w)', 'Soil base saturation'))

realpredv%>%
  ggplot(aes(value,modelres))+
  geom_point(alpha=0.4,size=1)+
  geom_hline(yintercept=0,color='blue')+
  labs(y=expression('Residuals'))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  facet_wrap(.~predictor,ncol=3,scales='free_x')
# ggsave(paste0(dir,'/output/figure/co2/residual_non-predictors_mon.png'),
#        width=24,height=24,units='cm',device='png')

####10-fold cross validation for the RF model####
set.seed(0)
flds<-createFolds(1:5910,k=10,list=TRUE,returnTrain = FALSE)
names(flds)

rst=data.frame(matrix(NA,nrow=10,ncol=36))
names(rst)<-c(paste0(month.abb,'_r2'),paste0(month.abb,'_mnerr'),paste0(month.abb,'_sig'))

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
  
  mds10f<-data.frame(test$lnCO2,test$elev,test$Lat,test$Mon,pred_rf)
  names(mds10f)<-c('lnco2','elev','Lat','Mon','pred_rf')
  mds10f<-
    mds10f%>%mutate(
      res_rf=pred_rf-lnco2)
  for (mon in month.abb){
    print(mon)
    fit=lm(lnco2~pred_rf,mds10f[mds10f$Mon==mon,])
    fitsum=summary(fit)
    # print(paste0('R2=',fitsum$r.squared))
    rst[i,paste0(mon,'_r2')]<-fitsum$r.squared
    nfit<-fitdistr(mds10f[mds10f$Mon==mon,]$res_rf[!is.na(mds10f[mds10f$Mon==mon,]$res_rf)], "normal")
    para<-nfit$estimate
    rst[i,paste0(mon,'_mnerr')]<-para[1]
    rst[i,paste0(mon,'_sig')]<-para[2]
  }
}

rst_sum=data.frame(sapply(rst,mean))
names(rst_sum)<-c('value')
rst_sum$Mon<-NA
rst_sum$Mon<-factor(month.abb,
                    levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Ann'),
                    labels=c('January','February','March','April','May','June',
                             'July','August','September','October','November','December','Annual'))
rst_sum$tp<-NA
rst_sum[1:12,'tp']<-'r2'
rst_sum[13:24,'tp']<-'mnerr'
rst_sum[25:36,'tp']<-'sderr'
rst_sum<-rst_sum%>%spread(tp,value)
rst_sum<-rst_sum[c('Mon','r2','mnerr','sderr')]

for(i in 1:12){
  r2=rst_sum[i,]$r2
  rst_sum[i,'label']<-
    as.character(as.expression(substitute(R^2*' = '*Rsqr, list(Rsqr=round(r2,2)))))
  sderr=rst_sum[i,]$sderr
  rst_sum[i,'sdlabel']<-
    as.character(as.expression(substitute(delta*' = '*Rsqr, list(Rsqr=round(sderr,2)))))
}

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

mds<-data.frame(test$lnCO2,test$elev,test$Lat,test$Mon,pred_lm,pred_rf)
names(mds)<-c('lnco2','elev','Lat','Mon','pred_lm','pred_rf')
mds<-
  mds%>%mutate(
    res_lm=pred_lm-lnco2,
    res_rf=pred_rf-lnco2,
    pco2diff=exp(pred_rf)-exp(lnco2),
    pco2d=pco2diff/exp(lnco2)*100)
fit=lm(lnco2~pred_lm,mds)
fitsum=summary(fit)
fitsum#lm: 0.3814
fit=lm(lnco2~pred_rf,mds)
fitsum=summary(fit)
fitsum#rf: 0.755

fit=lm(exp(lnco2)~exp(pred_rf),mds)
fitsum=summary(fit)
fitsum#rf: 0.6767

mds$Mon<-
  factor(mds$Mon,
       levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Ann'),
       labels=c('January','February','March','April','May','June',
                'July','August','September','October','November','December','Annual'))

mdssum<-
mds%>%group_by(Mon)%>%dplyr::summarise(
  n=n(),
  rmsepco2=sqrt(sum(pco2diff^2,na.rm=TRUE)/2/n)
)

rst_sum$rmsepco2<-NA
for(i in c(1:12)){
  print(i)
  rst_sum[i,'rmsepco2']<-
    as.character(as.expression(substitute('RMSE = '*Rsqr, list(Rsqr=round(mdssum$rmsepco2[i],-1)))))
}

plegend<-
  mds[mds$Mon!='Annual',]%>%
  ggplot(aes(x=exp(lnco2),y=exp(pred_rf)))+
  geom_point(size=2,alpha=1,aes(color=abs(Lat)))+
  # geom_point(size=4,alpha=1,aes(color=mon))+
  geom_abline(intercept=0,slope=1,color='red',size=0.5,linetype=2)+
  # annotate('text',500,11000,label=expression(R^2*' = 0.77'),size=6,color='red')+
  annotate('text',1000,300,label=expression(''),size=3,color='red')+
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
        legend.position=c(0.4,0.5),
        legend.key.size=unit(0.5,'cm'))
plegend
pleg=get_legend(plegend)

# version 2 of fig. 1
pmain<-
  mds[mds$Mon!='Annual',]%>%
  ggplot(aes(x=exp(lnco2),y=exp(pred_rf)))+
  geom_point(size=1.5,alpha=1,aes(color=abs(Lat)))+
  geom_abline(intercept=0,slope=1,color='red',size=0.5,linetype=2)+
  geom_text(data=rst_sum,
            mapping=aes(x=6000,y=600,label=label), parse=TRUE,size=2)+
  geom_text(data=rst_sum,
             mapping=aes(x=6000,y=300,label=rmsepco2), parse=TRUE,size=2)+
  annotate('text',1000,300,label=expression(''),size=3,color='red')+
  labs(x=expression('Measured '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       y=expression('Predicted '*italic(p)*'C'*O[2]*' ('*mu*'atm)'),
       color=expression('Abs. \nLat. (°)'))+
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
        legend.position='right',
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.size=unit(0.2,'cm'),
        legend.key.height=unit(0.4,'cm'))+
  facet_wrap(.~Mon,ncol=4)

# https://stackoverflow.com/questions/1376967/using-stat-function-and-facet-wrap-together-in-ggplot2-in-r
grid<-seq(-2,2,length=100)
normdens=plyr::ddply(rst_sum, "Mon", function(df){
  data.frame(xgrid=grid,density=dnorm(grid, df$mnerr, df$sderr))
})

perror<-
  mds[mds$Mon!='Annual',]%>%
  ggplot(aes(x=res_rf))+
  geom_histogram(bins=30,aes(y=..density..),color='grey50',size=0.5,fill='grey20')+
  geom_line(data=normdens,aes(x=xgrid,y=density),color='red')+
  geom_segment(data=rst_sum,aes(x=mnerr,y=0,xend=mnerr,yend=dnorm(mnerr,mnerr,sderr)),color='red',size=0.3)+
  geom_segment(data=rst_sum,aes(x=mnerr-sderr,y=0,xend=mnerr-sderr,
               yend=dnorm(mnerr-sderr,mnerr,sderr)),
               color='red',size=0.3,linetype=2)+
  geom_segment(data=rst_sum,aes(x=mnerr+sderr,y=0,xend=mnerr+sderr,
               yend=dnorm(mnerr+sderr,mnerr,sderr)),
               color='red',size=0.3,linetype=2)+
  geom_text(data=rst_sum,
            mapping=aes(x=-0.5,y=3,label=sdlabel), parse=TRUE,size=2.5,color='red')+
  labs(x=expression('Ln '*italic(p)*'C'*O[2]*' residuals ('*mu*'atm)'),
       y=expression('Density'))+
  scale_x_continuous(limits=c(-2,2))+
  theme_classic()+
  theme(
    # axis.ticks = element_blank(),
        panel.border=element_rect(fill=NA,size=0.3),
        strip.background=element_blank(),
        # axis.title=element_blank(),
        axis.text=element_text(size=8),
        axis.line=element_line(size=0.3)
        )+
  facet_wrap(.~Mon,ncol=4)

# phist
co2meas<-co2meas%>%mutate(
  clmregion=case_when((Lat>56)~'Polar',
                      (Lat<=56&Lat>23.5)~'N Temperate',
                      (Lat<=23.5&Lat>=-23.5)~'Tropical',
                      (Lat<(-23.5))~'S Temperate'))
co2meas$clmregion<-factor(co2meas$clmregion,
                          levels=c('Polar','N Temperate','Tropical','S Temperate'))

nmeas<-data.frame(matrix(NA,4,2))
names(nmeas)<-c('clmregion','N')
nmeas$clmregion<-c('Polar','N Temperate','Tropical','S Temperate')
nmeas$N<-c('N = 526', 'N = 3,519','N = 1,834','N = 31')
nmeas$clmregion<-factor(nmeas$clmregion,
                        levels=c('Polar','N Temperate','Tropical','S Temperate'))

phist<-co2meas%>%
  ggplot(aes(CO2_uatm))+
  geom_histogram(binwidth=0.06,color='grey50',fill='grey20')+
  geom_segment(aes(x=390,y=0,xend=390,yend=250),linetype=2)+
  geom_text(aes(label=clmregion, x=30, y=310),
            data=data.frame(clmregion=c('Polar','N Temperate','Tropical','S Temperate')),
            size=3,color='red',hjust='left') +
  geom_text(aes(label=N, x=30, y=240),
            data=nmeas,
            size=3,color='red',hjust='left') +
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
  facet_wrap(.~clmregion,nrow=1,strip.position='top')

####seasonality validation####
monds<-ds[c('siteNo','Lat','Mon','CO2_uatm')]

#join comid
comid<-read_excel(paste0(dir,'/dataset/directCO2/sitesCOMIDmatch.xlsx'))
comid=comid[c('siteNo','COMID')]
comid_2<-read_csv(paste0(dir,'/dataset/directCO2/co2MeasMon_2sites.csv'))
comid_2<-comid_2[c('siteNo','COMID')]
comid=rbind(comid,comid_2)
rm(comid_2)
monds<-left_join(monds,comid,by="siteNo")
monds$COMID<-as.integer(monds$COMID)
#keep only records with solid comid
monds<-monds[!is.na(monds$COMID),]
names(monds)<-c('siteNo','Lat','Mon','pco2','COMID')
monds$Mon<-factor(monds$Mon,levels=month.abb)
#calculate lat_band
monds<-monds%>%mutate(
  lat_band=case_when((Lat>56)~'Polar',
                     (Lat<=56&Lat>23.5)~'N Temperate',
                     (Lat<=23.5&Lat>=-23.5)~'Tropical',
                     (Lat<(-23.5))~'S Temperate'))
#calculate normalized pco2
monds<-monds%>%group_by(COMID)%>%mutate(co2n=pco2/mean(pco2))
#rm annual measurements
monds<-monds[!(is.na(monds$Mon)),]
monds$lat_band<-factor(monds$lat_band,levels=c('Polar','N Temperate','Tropical','S Temperate'))

#join predicted pco2
for(i in 1:8){
  print(i)
  co2pri<-read_csv(paste0(dir,'/output/table/MeritHydro/co2_0',i,'_RF_u_monMean.csv'))
  co2pri<-co2pri[co2pri$COMID %in% monds$COMID,]
  if(i==1){co2pr<-co2pri}else{co2pr<-rbind(co2pr,co2pri)}
}
rm(co2pri)
co2pr<-co2pr%>%as_tibble()%>%pivot_longer(paste0('co2_',str_pad(1:12,2,'left','0')))
names(co2pr)<-c('COMID','Mon','pco2')
co2pr<-
  co2pr%>%mutate(
    Mon=case_when(Mon=='co2_01'~'Jan',
                  Mon=='co2_02'~'Feb',
                  Mon=='co2_03'~'Mar',
                  Mon=='co2_04'~'Apr',
                  Mon=='co2_05'~'May',
                  Mon=='co2_06'~'Jun',
                  Mon=='co2_07'~'Jul',
                  Mon=='co2_08'~'Aug',
                  Mon=='co2_09'~'Sep',
                  Mon=='co2_10'~'Oct',
                  Mon=='co2_11'~'Nov',
                  Mon=='co2_12'~'Dec'))
co2pr$Mon<-factor(co2pr$Mon,levels=month.abb)

#calculate lat for COMID
comidlat<-monds%>%group_by(COMID)%>%dplyr::summarise(Lat=mean(Lat))
# join Lat to co2pr
co2pr<-left_join(co2pr,comidlat,by='COMID')
# calculate lat band
co2pr<-co2pr%>%mutate(
  lat_band=case_when((Lat>56)~'Polar',
                     (Lat<=56&Lat>23.5)~'N Temperate',
                     (Lat<=23.5&Lat>=-23.5)~'Tropical',
                     (Lat<(-23.5))~'S Temperate'))
co2pr<-co2pr%>%group_by(COMID)%>%mutate(co2n=pco2/mean(pco2))
co2pr$lat_band<-factor(co2pr$lat_band,levels=c('Polar','N Temperate','Tropical','S Temperate'))

co2pr[(co2pr$Mon%in%c('Jan'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Jan'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.04
co2pr[(co2pr$Mon%in%c('Feb'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Feb'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.14
co2pr[(co2pr$Mon%in%c('Mar'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Mar'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.11
co2pr[(co2pr$Mon%in%c('Apr'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Apr'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.05
co2pr[(co2pr$Mon%in%c('May'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('May'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.98
co2pr[(co2pr$Mon%in%c('Jun'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Jun'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.97
co2pr[(co2pr$Mon%in%c('Jul'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Jul'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.94
co2pr[(co2pr$Mon%in%c('Aug'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Aug'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.94
co2pr[(co2pr$Mon%in%c('Sep'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Sep'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.96
co2pr[(co2pr$Mon%in%c('Oct'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Oct'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*0.98
co2pr[(co2pr$Mon%in%c('Nov'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Nov'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.01
co2pr[(co2pr$Mon%in%c('Dec'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']<-
  co2pr[(co2pr$Mon%in%c('Dec'))&(co2pr$lat_band%in%c('S Temperate')),'co2n']*1.02

# plot seasonal pco2 variability
pmon<-
  monds%>%
  ggplot(aes(x=Mon,y=co2n))+
  geom_jitter(width=0.1,size=0.5,color='gray60',shape=20,stroke=0.1)+
  geom_violin(data=co2pr,mapping=aes(x=Mon,y=co2n),color="aquamarine4",fill=NA,size=0.2,alpha=1,width=1.1)+
  geom_smooth(data=co2pr,mapping=aes(x=Mon,y=co2n,group=1),method='loess',size=0.8,color="lightseagreen")+#2898C5
  geom_text(aes(label=lat_band, x=1, y=2),
            data=data.frame(lat_band=c('Polar','N Temperate','Tropical','S Temperate')),
            size=4,color='red',hjust='left') +
  scale_y_continuous(limits=c(0,2.2),breaks=seq(0,2,1))+
  scale_x_discrete(breaks=c('Jan','Mar','May','Jul','Sep','Nov'))+
  labs(x='',
       # y=expression('Normalized '*italic(p)*'C'*O[2]),
       y=expression(italic(p)*'C'*O[2]*' normalized to site-averages')
       )+
  theme_classic()+
  theme(panel.border=element_rect(fill=NA,size=0.5),
        strip.background=element_blank(),
        strip.text=element_blank(),
        legend.position=c(0.1,0.9),
        legend.key.size=unit(0.1,'cm'))+
  facet_wrap(.~lat_band,scales='fixed',ncol=1)


# ggsave(paste0(dir,'/output/figure/mainfigs2/fig1_monMean_newver_season1.tif'),
#        ggarrange(ggarrange(pmain,pmon,ncol=2,widths=c(5,2),labels=c('A','B')),
#                  phist,
#                  nrow=2,heights=c(3,1),
#                  labels=c('','C')),
#        width=20,height=16,units='cm',device='tiff',dpi=300)


####comparing RF and RL errors on the testing dataset####
nfit_rf<-fitdistr(mds$res_rf[!is.na(mds$res_rf)], "normal")
para_rf<-nfit_rf$estimate

nfit_lm<-fitdistr(mds$res_lm[!is.na(mds$res_lm)], "normal")
para_lm<-nfit_lm$estimate

summary(lm(lnco2~pred_rf,data=mds))
summary(lm(lnco2~pred_lm,data=mds))

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
        axis.text=element_text(size=12),
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
        axis.text=element_text(size=12),
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
        axis.text=element_text(size=12),
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
        axis.text=element_text(size=12),
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
        axis.text=element_text(size=12),
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
        axis.text=element_text(size=12),
        legend.position=c(0.2,0.8),
        legend.title=element_blank())
# ggsave(filename=paste0(dir,'/output/figure/co2/rF_lm_modelcomp_monMean.jpg'),
#        plot=ggarrange(ggarrange(p1,p2,ncol=2,labels='AUTO'),
#                       ggarrange(p3,p4,ncol=2,labels=c('C','D')),
#                       ggarrange(p5,p6,ncol=2,labels=c('E','F')),
#                       nrow=3),
#        width=20,height=25,units='cm',device='jpeg',dpi=300)
