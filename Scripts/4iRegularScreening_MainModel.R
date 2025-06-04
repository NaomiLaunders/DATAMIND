#Regular screening

####Any screening on same day####
rm(list = ls(all.names = TRUE))

####Libraries####
library(haven)
library(psych)
library(dplyr)
library(lubridate)
library(tidyverse)
library(tableone)
library(stringr)
library(tidyselect)
library(DescTools)
library(broom)
library(ivs)
library(sandwich)
library(lmtest)
library(nnet)
library(marginaleffects)
library(RStata)

####Load SMI cohort####

load("DatamindCPRD/Data/MainAnalysis.Rdata")
load("DatamindCPRD/Data/LipidAnalysis.Rdata")
load("DatamindCPRD/Data/BMIAnalysis.Rdata")
load("DatamindCPRD/Data/GlucAnalysis.Rdata")
load("DatamindCPRD/Data/BPAnalysis.Rdata")
load("DatamindCPRD/Data/AlcAnalysis.Rdata")
load("DatamindCPRD/Data/SmokeAnalysis.Rdata")

length(which(CPRD$start<=as.Date("2014-03-31") & CPRD$end==as.Date("2018-03-31")))
length(which(CPRD$start<=as.Date("2011-03-31") & CPRD$end>=as.Date("2014-03-31")))
length(which(CPRD$start<=as.Date("2004-03-31") & CPRD$end>=as.Date("2011-03-31")))

####Test Stata connection####
#chooseStataBin()
#Laptop
#options("RStata.StataPath" = "\"StataMP-64\"")
#options("RStata.StataVersion" = 18)

#Maple house
options("RStata.StataPath" = "\"StataMP-64\"")
options("RStata.StataVersion" = 17)

stata('di "Hello World"')
stata('net install regsave, from("https://raw.githubusercontent.com/reifjulian/regsave/master") replace')


####New variables####
CPRD$FU<-as.numeric((CPRD$end-CPRD$start)/365.25)
summary(CPRD$FU)

CPRD$Start04<-pmax(CPRD$start, as.Date("2004-04-01"))
CPRD$End04<-pmin(CPRD$end, as.Date("2011-03-31"))

CPRD$Start11<-pmax(CPRD$start, as.Date("2011-04-01"))
CPRD$End11<-pmin(CPRD$end, as.Date("2014-03-31"))

CPRD$Start14<-pmax(CPRD$start, as.Date("2014-04-01"))
CPRD$End14<-pmin(CPRD$end, as.Date("2018-03-31"))

CPRD$AgeAtStart04<-as.numeric(year(CPRD$Start04)-CPRD$yob)
CPRD$AgeAtStart11<-as.numeric(year(CPRD$Start11)-CPRD$yob)
CPRD$AgeAtStart14<-as.numeric(year(CPRD$Start14)-CPRD$yob)

CPRD$TimeSinceReg04<-CPRD$End04-CPRD$crd
CPRD$TimeSinceReg11<-CPRD$End11-CPRD$crd
CPRD$TimeSinceReg14<-CPRD$End14-CPRD$crd

CPRD$TimeSinceDiag04<-CPRD$End04-CPRD$FirstSMIDate
CPRD$TimeSinceDiag11<-CPRD$End11-CPRD$FirstSMIDate
CPRD$TimeSinceDiag14<-CPRD$End14-CPRD$FirstSMIDate

CPRD$TimeSinceReg04<-CPRD$End04-CPRD$crd
CPRD$TimeSinceReg11<-CPRD$End11-CPRD$crd
CPRD$TimeSinceReg14<-CPRD$End14-CPRD$crd

CPRD$FU04<-CPRD$End04-CPRD$Start04
CPRD$FU11<-CPRD$End11-CPRD$Start11
CPRD$FU14<-CPRD$End14-CPRD$Start14

CPRD$gender <- relevel(CPRD$gender, ref = "Female")
CPRD$LastDiag<-relevel(CPRD$LastDiag, ref = "bipolar")
CPRD$ethn_missing<-relevel(CPRD$ethn_missing, ref = "White")

#In period variables
#Exceptions to include non-QOF "physical health check exemption"
load("VariableExtracts/CPRD2018/MergedObs/SMIEx.Rdata")
PatExAll<-subset(PatExAll, patid %in% CPRD$patid & Group=="MH")

Dates<-select(CPRD, patid, start, end, Start04, End04, End11, Start11, End14, Start14, yob)
PatExAll<-merge(x=PatExAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatExAll<-subset(PatExAll, eventdate<=end & eventdate>yob)

Exempt2004<-subset(PatExAll, eventdate<=End04)
Exempt2011<-subset(PatExAll, eventdate<=End11)
Exempt2014<-subset(PatExAll, eventdate<=End14)

CPRD$EverMHEx2004<-0
CPRD$EverMHEx2004[CPRD$patid %in% Exempt2004$patid]<-1

CPRD$EverMHEx2011<-0
CPRD$EverMHEx2011[CPRD$patid %in% Exempt2011$patid]<-1

CPRD$EverMHEx2014<-0
CPRD$EverMHEx2014[CPRD$patid %in% Exempt2014$patid]<-1

#Other QOF
load("VariableExtracts/CPRD2018/MergedObs/OtherQOF.Rdata")
PatOtherQOFAll<-subset(PatOtherQOFAll, patid %in% CPRD$patid)

PatOtherQOFAll<-merge(x=PatOtherQOFAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatOtherQOFAll<-subset(PatOtherQOFAll, eventdate<=end & eventdate>yob)

OtherQOF2004<-subset(PatOtherQOFAll, eventdate<=End04)
OtherQOF2011<-subset(PatOtherQOFAll, eventdate<=End11)
OtherQOF2014<-subset(PatOtherQOFAll, eventdate<=End14)

CPRD$EverOtherQOF2004<-0
CPRD$EverOtherQOF2004[CPRD$patid %in% OtherQOF2004$patid]<-1

CPRD$EverOtherQOF2011<-0
CPRD$EverOtherQOF2011[CPRD$patid %in% OtherQOF2011$patid]<-1

CPRD$EverOtherQOF2014<-0
CPRD$EverOtherQOF2014[CPRD$patid %in% OtherQOF2014$patid]<-1

#Antipsych
load("VariableExtracts/CPRD2018/MergedObs/APProd.Rdata")
PatAPAll<-subset(PatAPAll, patid %in% CPRD$patid)
PatAPAll<-merge(x=PatAPAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatAPAll<-subset(PatAPAll, eventdate>=start & eventdate<=end)

AP2004<-subset(PatAPAll, eventdate<=End04)
AP2011<-subset(PatAPAll, eventdate<=End11)
AP2014<-subset(PatAPAll, eventdate<=End14)

CPRD$AP2004<-0
CPRD$AP2004[CPRD$patid %in% AP2004$patid]<-1

CPRD$AP2011<-0
CPRD$AP2011[CPRD$patid %in% AP2011$patid]<-1

CPRD$AP2014<-0
CPRD$AP2014[CPRD$patid %in% AP2014$patid]<-1

####Limit look ups####

Alc<-subset(Alc, patid %in% CPRD$patid)
BMI<-subset(BMI, patid %in% CPRD$patid)
Smoke<-subset(Smoke, patid %in% CPRD$patid)
Lipid<-subset(Lipid, patid %in% CPRD$patid)
Gluc<-subset(Gluc, patid %in% CPRD$patid)
BP<-subset(BP, patid %in% CPRD$patid)

####Full cohort####
Active<-select(CPRD, patid, Act0405, Act0506, Act0607, Act0708, Act0809, Act0910, Act1011, Act1112, Act1213, Act1314, Act1415, Act1516, Act1617, Act1718)
#Alcohol
Alc<-merge(x=Alc, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
Alc[,6:19]<-lapply(Alc[,6:19], as.character)
Alc[,6:19]<-lapply(Alc[,6:19], as.numeric)  


AlcYearAll<-Alc%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(AlcYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(AlcYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  AlcYearAll[paste0(ScreenVar[i])]<-if_else(AlcYearAll[paste0(ActVar[i])]==0, 0, if_else(AlcYearAll[paste0(ActVar[i])]==1 & AlcYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  AlcYearAll[paste0(ScreenVar[i])][is.na(AlcYearAll[paste0(ScreenVar[i])])]<-0
}

AlcYearAll<-AlcYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(AlcComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         AlcComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         AlcComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))
  
AlcYearAll<-select(AlcYearAll, patid, AlcComplete04, AlcComplete11, AlcComplete14)

#BMI
BMI<-merge(x=BMI, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
BMI[,6:19]<-lapply(BMI[,6:19], as.character)
BMI[,6:19]<-lapply(BMI[,6:19], as.numeric)  

BMIYearAll<-BMI%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(BMIYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(BMIYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  BMIYearAll[paste0(ScreenVar[i])]<-if_else(BMIYearAll[paste0(ActVar[i])]==0, 0, if_else(BMIYearAll[paste0(ActVar[i])]==1 & BMIYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  BMIYearAll[paste0(ScreenVar[i])][is.na(BMIYearAll[paste0(ScreenVar[i])])]<-0
}

BMIYearAll<-BMIYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(BMIComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         BMIComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         BMIComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))

BMIYearAll<-select(BMIYearAll, patid, BMIComplete04, BMIComplete11, BMIComplete14)

#Smoking
Smoke<-merge(x=Smoke, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
Smoke[,6:19]<-lapply(Smoke[,6:19], as.character)
Smoke[,6:19]<-lapply(Smoke[,6:19], as.numeric)  

SmokeYearAll<-Smoke%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(SmokeYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(SmokeYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  SmokeYearAll[paste0(ScreenVar[i])]<-if_else(SmokeYearAll[paste0(ActVar[i])]==0, 0, if_else(SmokeYearAll[paste0(ActVar[i])]==1 & SmokeYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  SmokeYearAll[paste0(ScreenVar[i])][is.na(SmokeYearAll[paste0(ScreenVar[i])])]<-0
}

SmokeYearAll<-SmokeYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(SmokeComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         SmokeComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         SmokeComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))

SmokeYearAll<-select(SmokeYearAll, patid, SmokeComplete04, SmokeComplete11, SmokeComplete14)

#Lipids
Lipid<-merge(x=Lipid, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
Lipid[,6:19]<-lapply(Lipid[,6:19], as.character)
Lipid[,6:19]<-lapply(Lipid[,6:19], as.numeric)  

LipidYearAll<-Lipid%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(LipidYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(LipidYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  LipidYearAll[paste0(ScreenVar[i])]<-if_else(LipidYearAll[paste0(ActVar[i])]==0, 0, if_else(LipidYearAll[paste0(ActVar[i])]==1 & LipidYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  LipidYearAll[paste0(ScreenVar[i])][is.na(LipidYearAll[paste0(ScreenVar[i])])]<-0
}

LipidYearAll<-LipidYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(LipidComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         LipidComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         LipidComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))

LipidYearAll<-select(LipidYearAll, patid, LipidComplete04, LipidComplete11, LipidComplete14)

#Glucose
Gluc<-merge(x=Gluc, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
Gluc[,6:19]<-lapply(Gluc[,6:19], as.character)
Gluc[,6:19]<-lapply(Gluc[,6:19], as.numeric)  

GlucYearAll<-Gluc%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(GlucYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(GlucYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  GlucYearAll[paste0(ScreenVar[i])]<-if_else(GlucYearAll[paste0(ActVar[i])]==0, 0, if_else(GlucYearAll[paste0(ActVar[i])]==1 & GlucYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  GlucYearAll[paste0(ScreenVar[i])][is.na(GlucYearAll[paste0(ScreenVar[i])])]<-0
}

GlucYearAll<-GlucYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(GlucComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         GlucComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         GlucComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))

GlucYearAll<-select(GlucYearAll, patid, GlucComplete04, GlucComplete11, GlucComplete14)

#BP
BP<-merge(x=BP, y=Active, by="patid", all.x=TRUE, all.y=FALSE)
BP[,6:19]<-lapply(BP[,6:19], as.character)
BP[,6:19]<-lapply(BP[,6:19], as.numeric)  

BPYearAll<-BP%>%
  mutate(Year=case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "y0405",
                        eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "y0506",
                        eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "y0607",
                        eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "y0708",
                        eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "y0809",
                        eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "y0910",
                        eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "y1011",
                        eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "y1112",
                        eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "y1213",
                        eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "y1314",
                        eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "y1415",
                        eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "y1516",
                        eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "y1617",
                        eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "y1718",
                        TRUE ~ "HELP"))%>%
  subset(Year!="HELP")%>%
  select(-eventdate)%>%
  distinct()%>%
  group_by(patid)%>%
  mutate(n=1)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill = 0)

#Set as not screened if not active during that time
ActVar<-vars_select(names(BPYearAll), starts_with('Act', ignore.case = TRUE))
ScreenVar<-vars_select(names(BPYearAll), starts_with('y', ignore.case = TRUE))

for (i in (1:length(ScreenVar))) {
  BPYearAll[paste0(ScreenVar[i])]<-if_else(BPYearAll[paste0(ActVar[i])]==0, 0, if_else(BPYearAll[paste0(ActVar[i])]==1 & BPYearAll[paste0(ScreenVar[i])]>0, 1, 0))
  BPYearAll[paste0(ScreenVar[i])][is.na(BPYearAll[paste0(ScreenVar[i])])]<-0
}

BPYearAll<-BPYearAll%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(5:11)),#Count up the number active years
         Active11=sum(c_across(12:14)),
         Active14=sum(c_across(15:18)),
         Screened04=sum(c_across(19:25)),
         Screened11=sum(c_across(26:28)),
         Screened14=sum(c_across(29:32)))%>%
  mutate(BPComplete04=case_when(Screened04==0 ~ "Never",
                                 Active04==Screened04 ~ "Complete",
                                 Active04>Screened04 ~ "Mixed",
                                 TRUE ~ "Check"),
         BPComplete11=case_when(Screened11==0 ~ "Never",
                                 Active11==Screened11 ~ "Complete",
                                 Active11>Screened11 ~ "Mixed",
                                 TRUE ~ "Check"),
         BPComplete14=case_when(Screened14==0 ~ "Never",
                                 Active14==Screened14 ~ "Complete",
                                 Active14>Screened14 ~ "Mixed",
                                 TRUE ~ "Check"))

BPYearAll<-select(BPYearAll, patid, BPComplete04, BPComplete11, BPComplete14)

#Merge back to data
CPRD<-merge(x=CPRD, y=AlcYearAll, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-merge(x=CPRD, y=BMIYearAll, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-merge(x=CPRD, y=BPYearAll, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-merge(x=CPRD, y=LipidYearAll, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-merge(x=CPRD, y=GlucYearAll, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-merge(x=CPRD, y=SmokeYearAll, by="patid", all.x=TRUE, all.y=FALSE)

####Set active as numeric, and blank Alc to smoke as never
which(colnames(CPRD) == "Act0001")
which(colnames(CPRD) == "Act1718")

CPRD[,125:142]<-lapply(CPRD[,125:142], as.character)
CPRD[,125:142]<-lapply(CPRD[,125:142], as.numeric)  

which(colnames(CPRD) == "AlcComplete04")
which(colnames(CPRD) == "SmokeComplete14")

CPRD <- CPRD %>% 
  mutate_at(c(204:221), ~replace_na(.,"Never"))

####Split those who have at least 90 days active in each year
which(colnames(CPRD) == "Act0405")
which(colnames(CPRD) == "Act1011")
which(colnames(CPRD) == "Act1314")
which(colnames(CPRD) == "Act1718")

CPRDMod<-CPRD%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Act04 = sum(c_across(129:135)),
         Act11 = sum(c_across(136:138)),
         Act18 = sum(c_across(139:142)))

length(which(CPRDMod$Act04>=1 | CPRDMod$Act11>=1|CPRDMod$Act18>=1))

CPRD2004<-subset(CPRDMod, Act04>=1)
CPRD2011<-subset(CPRDMod, Act11>=1)
CPRD2014<-subset(CPRDMod, Act18>=1)

CPRDCheck<-subset(CPRD, patid %in% CPRD2004$patid |patid %in% CPRD2011$patid |patid %in% CPRD2014$patid)

CPRD2004<-CPRD2004%>%
  mutate(AllComplete04=case_when(SmokeComplete04=="Complete" & AlcComplete04=="Complete"& BMIComplete04=="Complete"
                                 & LipidComplete04=="Complete" & GlucComplete04=="Complete"& BPComplete04=="Complete" ~ "Complete",
                                 SmokeComplete04=="Never" & AlcComplete04=="Never"& BMIComplete04=="Never"
                                 & LipidComplete04=="Never" & GlucComplete04=="Never"& BPComplete04=="Never" ~ "Never",
                                 TRUE ~ "Mixed")) 
CPRD2011<-CPRD2011%>%
  mutate(AllComplete11=case_when(SmokeComplete11=="Complete" & AlcComplete11=="Complete"& BMIComplete11=="Complete"
                                 & LipidComplete11=="Complete" & GlucComplete11=="Complete"& BPComplete11=="Complete" ~ "Complete",
                                 SmokeComplete11=="Never" & AlcComplete11=="Never"& BMIComplete11=="Never"
                                 & LipidComplete11=="Never" & GlucComplete11=="Never"& BPComplete11=="Never" ~ "Never",
                                 TRUE ~ "Mixed"))
CPRD2014<-CPRD2014%>%
  mutate(AllComplete14=case_when(SmokeComplete14=="Complete" & AlcComplete14=="Complete"& BMIComplete14=="Complete"
                                 & LipidComplete14=="Complete" & GlucComplete14=="Complete"& BPComplete14=="Complete" ~ "Complete",
                                 SmokeComplete14=="Never" & AlcComplete14=="Never"& BMIComplete14=="Never"
                                 & LipidComplete14=="Never" & GlucComplete14=="Never"& BPComplete14=="Never" ~ "Never",
                                 TRUE ~ "Mixed"))



CPRD2004[,204:221]<-lapply(CPRD2004[,204:221], factor, levels=c("Mixed", "Complete", "Never"))
CPRD2011[,204:221]<-lapply(CPRD2011[,204:221], factor, levels=c("Mixed", "Complete", "Never"))
CPRD2014[,204:221]<-lapply(CPRD2014[,204:221], factor, levels=c("Mixed", "Complete", "Never"))

CPRD2004[,225]<-lapply(CPRD2004[,225], factor, levels=c("Mixed", "Complete", "Never"))
CPRD2011[,225]<-lapply(CPRD2011[,225], factor, levels=c("Mixed", "Complete", "Never"))
CPRD2014[,225]<-lapply(CPRD2014[,225], factor, levels=c("Mixed", "Complete", "Never"))

#year of end of follow up
CPRD2004<-CPRD2004%>%
  mutate(Year04 = case_when(End04>="2004-04-01" & End04<="2005-03-31" ~ "0405",
                            End04>="2005-04-01" & End04<="2006-03-31" ~ "0506",
                            End04>="2006-04-01" & End04<="2007-03-31" ~ "0607",
                            End04>="2007-04-01" & End04<="2008-03-31" ~ "0708",
                            End04>="2008-04-01" & End04<="2009-03-31" ~ "0809",
                            End04>="2009-04-01" & End04<="2010-03-31" ~ "0910",
                            End04>="2010-04-01" & End04<="2011-03-31" ~ "1011",
                               TRUE ~ "HELP"))
table(CPRD2004$Year04)

CPRD2011<-CPRD2011%>%
  mutate(Year11 = case_when(End11>="2011-04-01" & End11<="2012-03-31" ~ "1112",
                            End11>="2012-04-01" & End11<="2013-03-31" ~ "1213",
                            End11>="2013-04-01" & End11<="2014-03-31" ~ "1314",
                            TRUE ~ "HELP"))
table(CPRD2011$Year11)

CPRD2014<-CPRD2014%>%
  mutate(Year14 = case_when(End14>="2014-04-01" & End14<="2015-03-31" ~ "1415",
                            End14>="2015-04-01" & End14<="2016-03-31" ~ "1516",
                            End14>="2016-04-01" & End14<="2017-03-31" ~ "1617",
                            End14>="2017-04-01" & End14<="2018-03-31" ~ "1718",
                            TRUE ~ "HELP"))
table(CPRD2014$Year14)

CPRD2004$Year04<-as.factor(CPRD2004$Year04)
CPRD2011$Year11<-as.factor(CPRD2011$Year11)
CPRD2014$Year14<-as.factor(CPRD2014$Year14)

CPRD2004$Year04 <- relevel(CPRD2004$Year04, ref = "1011")
CPRD2011$Year11 <- relevel(CPRD2011$Year11, ref = "1314")
CPRD2014$Year14 <- relevel(CPRD2014$Year14, ref = "1718")

CPRD2004$FullIMD <- relevel(CPRD2004$FullIMD, ref = "1")
CPRD2011$FullIMD <- relevel(CPRD2011$FullIMD, ref = "1")
CPRD2014$FullIMD <- relevel(CPRD2014$FullIMD, ref = "1")

CPRD2004[,204:221]<-lapply(CPRD2004[,204:221], relevel, ref="Mixed")
CPRD2011[,204:221]<-lapply(CPRD2011[,204:221], relevel, ref="Mixed")
CPRD2014[,204:221]<-lapply(CPRD2014[,204:221], relevel, ref="Mixed")

CPRD2004[,225]<-lapply(CPRD2004[,225], relevel, ref="Mixed")
CPRD2011[,225]<-lapply(CPRD2011[,225], relevel, ref="Mixed")
CPRD2014[,225]<-lapply(CPRD2014[,225], relevel, ref="Mixed")

#Currently have one year follow up and are active in that year if active for last 90 days
#Patients in only one year can only be complete or never screened so need them to be active in at least 2 years

CPRD2004$Incl<-0
CPRD2004$Incl[CPRD2004$Act04>=2]<-1

CPRD2011$Incl<-0
CPRD2011$Incl[CPRD2011$Act11>=2]<-1

CPRD2014$Incl<-0
CPRD2014$Incl[CPRD2014$Act18>=2]<-1

#Check 2004
MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Incl", data=CPRD2004, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Full2004.csv")

#ANd for 2011

MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Incl", data=CPRD2011, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Full2011.csv")

#ANd for 2014

MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Incl", data=CPRD2014, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Full2014.csv")

CPRD2004_2<-subset(CPRD2004, Incl==1)
CPRD2011_2<-subset(CPRD2011, Incl==1)
CPRD2014_2<-subset(CPRD2014, Incl==1)

####How many patients are in all groups?####
CPRDCheck<-subset(CPRD, patid %in% CPRD2004_2$patid | patid %in% CPRD2011_2$patid |patid %in% CPRD2014_2$patid)
CPRDCheck1<-subset(CPRD, patid %in% CPRD2011_2$patid |patid %in% CPRD2014_2$patid)
###Multinomial logistic regression####

####Unadjusted####

#Age

stata('mlogit AllComplete04 AgeAtStart04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Age04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 AgeAtStart11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Age11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 AgeAtStart14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Age14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Sex
stata('mlogit AllComplete04  i.gender, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Sex04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.gender, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Sex11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.gender, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Sex14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Ethnicity
stata('mlogit AllComplete04  i.ethn_missing , baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Ethn04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.ethn_missing, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Ethn11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.ethn_missing, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Ethn14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Country
stata('mlogit AllComplete04  i.country , baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Country04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.country, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Country11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.country, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Country14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#SMI diagnosis
stata('mlogit AllComplete04  i.LastDiag, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/SMI04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.LastDiag, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/SMI11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.LastDiag, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/SMI14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Exceptions
stata('mlogit AllComplete04  i.EverMHEx2004, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Ex04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.EverMHEx2011, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Ex11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.EverMHEx2014, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Ex14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#QOF
stata('mlogit AllComplete04  i.EverOtherQOF2004 , baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/QOF04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.EverOtherQOF2011, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/QOF11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.EverOtherQOF2014, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/QOF14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#APs
stata('mlogit AllComplete04  i.AP2004, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/AP04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 i.AP2011, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/AP11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 i.AP2014, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/AP14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Since diag
stata('mlogit AllComplete04 TimeSinceDiag04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/TimeDiag04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 TimeSinceDiag11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/TimeDiag11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 TimeSinceDiag14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/TimeDiag14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Since reg
stata('mlogit AllComplete04 TimeSinceReg04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/TimeReg04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 TimeSinceReg11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/TimeReg11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 TimeSinceReg14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/TimeReg14.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Since FU
stata('mlogit AllComplete04  FU04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/FU04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11  FU11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/FU11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14  FU14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/FU14.csv", replace ci nostar eform', data.in = CPRD2014_2)


#Year of record
stata('mlogit AllComplete04  i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/Year04.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11  i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/Year11.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14  i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/Year14.csv", replace ci nostar eform', data.in = CPRD2014_2)


####Mutually adjusted####



#All
stata('mlogit AllComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/All04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AllComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/All11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AllComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/All14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)


#Smoking
stata('mlogit SmokeComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store smoke04
      esttab smoke04 using "DatamindCPRD/Outputs/smoke04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit SmokeComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store smoke11
  esttab smoke11 using "DatamindCPRD/Outputs/smoke11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit SmokeComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store smoke14
  esttab smoke14 using "DatamindCPRD/Outputs/smoke14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)


#Alcohol
stata('mlogit AlcComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Alc04
      esttab Alc04 using "DatamindCPRD/Outputs/Alc04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit AlcComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Alc11
  esttab Alc11 using "DatamindCPRD/Outputs/Alc11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit AlcComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Alc14
  esttab Alc14 using "DatamindCPRD/Outputs/Alc14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)

#BMI
stata('mlogit BMIComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store BMI04
      esttab BMI04 using "DatamindCPRD/Outputs/BMI04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit BMIComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BMI11
  esttab BMI11 using "DatamindCPRD/Outputs/BMI11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit BMIComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BMI14
  esttab BMI14 using "DatamindCPRD/Outputs/BMI14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Lipids
stata('mlogit LipidComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Lipids04
      esttab Lipids04 using "DatamindCPRD/Outputs/Lipids04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit LipidComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Lipids11
  esttab Lipids11 using "DatamindCPRD/Outputs/Lipids11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit LipidComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Lipids14
  esttab Lipids14 using "DatamindCPRD/Outputs/Lipids14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)

#Glucose
stata('mlogit GlucComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Gluc04
      esttab Gluc04 using "DatamindCPRD/Outputs/Gluc04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit GlucComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Gluc11
  esttab Gluc11 using "DatamindCPRD/Outputs/Gluc11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit GlucComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Gluc14
  esttab Gluc14 using "DatamindCPRD/Outputs/Gluc14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)

#BP
stata('mlogit BPComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store BP04
      esttab BP04 using "DatamindCPRD/Outputs/BP04Full.csv", replace ci nostar eform', data.in = CPRD2004_2)

stata('mlogit BPComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BP11
  esttab BP11 using "DatamindCPRD/Outputs/BP11Full.csv", replace ci nostar eform', data.in = CPRD2011_2)

stata('mlogit BPComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BP14
  esttab BP14 using "DatamindCPRD/Outputs/BP14Full.csv", replace ci nostar eform', data.in = CPRD2014_2)

####IMD sensitivity####
#All
IMD2004<-subset(CPRD2004_2, !is.na(FullIMD))
IMD2011<-subset(CPRD2011_2, !is.na(FullIMD))
IMD2014<-subset(CPRD2014_2, !is.na(FullIMD))

####How many patients are in all groups?####
IMDCheck<-subset(CPRD, patid %in% IMD2004$patid | patid %in% IMD2011$patid |patid %in% IMD2014$patid)
IMDCheck1<-subset(CPRD, patid %in% IMD2011$patid |patid %in% IMD2014$patid)

#All
stata('mlogit AllComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/All04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit AllComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/All11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit AllComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/All14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

#Smoking
stata('mlogit SmokeComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store smoke04
      esttab smoke04 using "DatamindCPRD/Outputs/smoke04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit SmokeComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store smoke11
  esttab smoke11 using "DatamindCPRD/Outputs/smoke11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit SmokeComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store smoke14
  esttab smoke14 using "DatamindCPRD/Outputs/smoke14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)


#Alcohol
stata('mlogit AlcComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Alc04
      esttab Alc04 using "DatamindCPRD/Outputs/Alc04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit AlcComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Alc11
  esttab Alc11 using "DatamindCPRD/Outputs/Alc11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit AlcComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Alc14
  esttab Alc14 using "DatamindCPRD/Outputs/Alc14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

#BMI
stata('mlogit BMIComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store BMI04
      esttab BMI04 using "DatamindCPRD/Outputs/BMI04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit BMIComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BMI11
  esttab BMI11 using "DatamindCPRD/Outputs/BMI11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit BMIComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BMI14
  esttab BMI14 using "DatamindCPRD/Outputs/BMI14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

#Lipids
stata('mlogit LipidComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Lipids04
      esttab Lipids04 using "DatamindCPRD/Outputs/Lipids04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit LipidComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Lipids11
  esttab Lipids11 using "DatamindCPRD/Outputs/Lipids11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit LipidComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Lipids14
  esttab Lipids14 using "DatamindCPRD/Outputs/Lipids14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

#Glucose
stata('mlogit GlucComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store Gluc04
      esttab Gluc04 using "DatamindCPRD/Outputs/Gluc04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit GlucComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Gluc11
  esttab Gluc11 using "DatamindCPRD/Outputs/Gluc11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit GlucComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store Gluc14
  esttab Gluc14 using "DatamindCPRD/Outputs/Gluc14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

#BP
stata('mlogit BPComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04 FU04 i.Year04 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
      estimates store BP04
      esttab BP04 using "DatamindCPRD/Outputs/BP04IMDFull.csv", replace ci nostar eform', data.in = IMD2004)

stata('mlogit BPComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11 FU11 i.Year11 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BP11
  esttab BP11 using "DatamindCPRD/Outputs/BP11IMDFull.csv", replace ci nostar eform', data.in = IMD2011)

stata('mlogit BPComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14 FU14 i.Year14 i.FullIMD, baseoutcome(1) vce(cluster pracid) rrr
  estimates store BP14
  esttab BP14 using "DatamindCPRD/Outputs/BP14IMDFull.csv", replace ci nostar eform', data.in = IMD2014)

####Complete case sensitivity####
CPRD2004$Complete<-0
CPRD2004$Complete[CPRD2004$Act04==7]<-1

CPRD2011$Complete<-0
CPRD2011$Complete[CPRD2011$Act11==3]<-1

CPRD2014$Complete<-0
CPRD2014$Complete[CPRD2014$Act18==4]<-1

#Check 2004
MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Complete", data=CPRD2004, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Comp2004.csv")

#And for 2011

MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Complete", data=CPRD2011, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Comp2011.csv")

#ANd for 2014

MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "died", "AgeAtDeath",  "LastDiag", "EverMHEx", "AnyPres")

Table1<-CreateTableOne(vars=MyVars,  strata="Complete", data=CPRD2014, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Sup_Comp2014.csv")

CPRD2004Comp<-subset(CPRD2004, Complete==1)
CPRD2011Comp<-subset(CPRD2011, Complete==1)
CPRD2014Comp<-subset(CPRD2014, Complete==1)

####How many patients are in all groups?####
CompCheck<-subset(CPRD, patid %in% CPRD2004Comp$patid | patid %in% CPRD2011Comp$patid |patid %in% CPRD2014Comp$patid)
CompCheck1<-subset(CPRD, patid %in% CPRD2011Comp$patid |patid %in% CPRD2014Comp$patid)
stata('mlogit AllComplete04 AgeAtStart04 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2004 i.EverOtherQOF2004 i.AP2004 TimeSinceDiag04 TimeSinceReg04, baseoutcome(1) vce(cluster pracid) rrr
      estimates store All04
      esttab All04 using "DatamindCPRD/Outputs/All04Comp.csv", replace ci nostar eform', data.in = CPRD2004Comp)

stata('mlogit AllComplete11 AgeAtStart11 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2011 i.EverOtherQOF2011 i.AP2011 TimeSinceDiag11 TimeSinceReg11, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All11
  esttab All11 using "DatamindCPRD/Outputs/All11Comp.csv", replace ci nostar eform', data.in = CPRD2011Comp)

stata('mlogit AllComplete14 AgeAtStart14 i.gender i.ethn_missing i.country i.LastDiag i.EverMHEx2014 i.EverOtherQOF2014 i.AP2014 TimeSinceDiag14 TimeSinceReg14, baseoutcome(1) vce(cluster pracid) rrr
  estimates store All14
  esttab All14 using "DatamindCPRD/Outputs/All14Comp.csv", replace ci nostar eform', data.in = CPRD2014Comp)



#How many are 
MyVars<-c("AllComplete04", "AlcComplete04", "BPComplete04", "SmokeComplete04", "LipidComplete04", "GlucComplete04", "BMIComplete04")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD2004_2, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Count_Full04.csv")

MyVars<-c("AllComplete11", "AlcComplete11", "BPComplete11", "SmokeComplete11", "LipidComplete11", "GlucComplete11", "BMIComplete11")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD2011_2, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Count_Full11.csv")

MyVars<-c("AllComplete14", "AlcComplete14", "BPComplete14", "SmokeComplete14", "LipidComplete14", "GlucComplete14", "BMIComplete14")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD2014_2, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Count_Full14.csv")

table(CPRD2004Comp$AllComplete04)
table(CPRD2011Comp$AllComplete11)
table(CPRD2014Comp$AllComplete14)

prop.table(table(CPRD2011Comp$AllComplete11))*100
prop.table(table(CPRD2014Comp$AllComplete14))*100

table(IMD2004$AllComplete04)
table(IMD2011$AllComplete11)
table(IMD2014$AllComplete14)



