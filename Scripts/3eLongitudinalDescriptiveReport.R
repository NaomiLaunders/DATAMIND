rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(lubridate)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/FinalCohortActiveFactors.Rdata")

####BMI####
load("VariableExtracts/CPRD2018/MergedObs/BMIValue.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/BMICalc.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/BMIAll.Rdata")

BMICalc<-select(BMICalc, patid, eventdate)
BMIValue<-select(BMIValue, patid, eventdate)
BMICat<-subset(PatBMIAll, !is.na(BMICat))
BMICat<-select(BMICat, patid, eventdate)

BMI<-rbind(BMICalc, BMIValue, BMICat)

rm(BMICalc, BMIValue, BMICat, PatBMIAll)

BMI<-subset(BMI, patid %in% CPRD$patid)
length(unique(BMI$patid))

BMI<-subset(BMI, eventdate>="2000-04-01" & eventdate<="2018-03-31")
length(unique(BMI$patid))

Dates<-select(CPRD, patid, start, end, country)
BMI<-merge(x=BMI, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMI<-subset(BMI, eventdate>=start & eventdate<=end)
length(unique(BMI$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-BMI%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

rm(BMI)
#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

167641/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med)

####Total chol####

load("VariableExtracts/CPRD2018/MergedObs/AllLipids_FinalinclCat.Rdata")

AllLipids<-subset(AllLipids, patid %in% CPRD$patid)

length(unique(AllLipids$patid))
Lipids<-subset(AllLipids, !is.na(tchol) | !is.na(Cat_totalcholesterol))
length(unique(Lipids$patid))

Lipids<-subset(Lipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Lipids<-merge(x=Lipids, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipids<-subset(Lipids, eventdate>=start & eventdate<=end)
length(unique(Lipids$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Lipids%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

135736/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Scotland
Scot<-Lipids%>%
  subset(country=="Scotland")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/14010*100)

#Wales
Wales<-Lipids%>%
  subset(country=="Wales")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/11841*100)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med)

####TC:HDL####
Lipids<-subset(AllLipids, !is.na(tchdlratio)| !is.na(Cat_tchdlratio))
length(unique(Lipids$patid))

Lipids<-subset(Lipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Lipids<-merge(x=Lipids, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipids<-subset(Lipids, eventdate>=start & eventdate<=end)
length(unique(Lipids$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Lipids%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

126162/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Scotland
Scot<-Lipids%>%
  subset(country=="Scotland")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/14010*100)

#Wales
Wales<-Lipids%>%
  subset(country=="Wales")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/11841*100)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, Lipids, AllLipids)

####Blood pressure####

load("VariableExtracts/CPRD2018/MergedObs/AllBP_FinalinclCat.Rdata")

AllBP<-subset(AllBP, patid %in% CPRD$patid)

length(unique(AllBP$patid))

BP<-subset(AllBP, eventdate>="2000-04-01" & eventdate<="2018-03-31")

BP<-merge(x=BP, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BP<-subset(BP, eventdate>=start & eventdate<=end)
length(unique(BP$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-BP%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

189557/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, AllBP)

####Require both BP####

LimitBP<-subset(BP, (!is.na(Diastolic)&!is.na(Systolic)) | (!is.na(Cat_Diastolic)&!is.na(Cat_Systolic)))
length(unique(LimitBP$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-LimitBP%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

189162/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Scotland
Scot<-LimitBP%>%
  subset(country=="Scotland")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/14010*100)

#Wales
Wales<-LimitBP%>%
  subset(country=="Wales")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/11841*100)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, BP, LimitBP)

####Blood glucose####

load("VariableExtracts/CPRD2018/MergedObs/AllGlucose_FinalinclCat.Rdata")

AllGlucose<-subset(AllGlucose, patid %in% CPRD$patid)

length(unique(AllGlucose$patid))

Gluc<-subset(AllGlucose, !is.na(Glucose) | !is.na(Cat_Glucose))

Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)
length(unique(Gluc$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Gluc%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

136374/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Scotland
Scot<-Gluc%>%
  subset(country=="Scotland")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/14010*100)

#Wales
Wales<-Gluc%>%
  subset(country=="Wales")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/11841*100)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, Gluc)

####HbA1c####
Gluc<-subset(AllGlucose, !is.na(Hba1c) | !is.na(Cat_Hba1c))
Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)
length(unique(Gluc$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Gluc%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

86556/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Scotland
Scot<-Gluc%>%
  subset(country=="Scotland")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/14010*100)

#Wales
Wales<-Gluc%>%
  subset(country=="Wales")%>%
  group_by(patid)%>%
  summarise(n=n())%>%
  ungroup()%>%
  summarise(n=n(), Perc=n/11841*100)
#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, Gluc, AllGlucose)
####Smoking####

load("VariableExtracts/CPRD2018/MergedObs/Smoke_FinalCat.Rdata")

CatSmoke<-subset(CatSmoke, patid %in% CPRD$patid)

length(unique(CatSmoke$patid))

Smoke<-subset(CatSmoke, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Smoke<-merge(x=Smoke, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Smoke<-subset(Smoke, eventdate>=start & eventdate<=end)
length(unique(Smoke$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Smoke%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto6mo=case_when(eventdate-First>7 & eventdate-First<=182.625 ~ 1, TRUE ~ 0),
         upto1yr=case_when(eventdate-First>7 & eventdate-First<=365.25 ~ 1, TRUE ~ 0),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0),
         ever=case_when(eventdate-First>0 ~ 1, TRUE ~ 0 ))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto6mo=sum(upto6mo), upto1yr=sum(upto1yr), upto2yr=sum(upto2yr), ever=sum(ever))%>%
  mutate(upto6mo=case_when(upto6mo>0~1, TRUE ~0),
         upto1yr=case_when(upto1yr>0~1, TRUE ~0),
         upto2yr=case_when(upto2yr>0~1, TRUE ~0),
         ever=case_when(ever>0~1, TRUE ~0))

190933/216136

table(LongSum$upto6mo)
table(LongSum$upto1yr)
table(LongSum$upto2yr)
table(LongSum$ever)

#Medians

Med<-Long%>%
  ungroup()%>%
  mutate(Length=eventdate-First)%>%
  subset(Length>0)%>%
  group_by(patid)%>%
  mutate(Next=min(Length))%>%
  select(patid, Next)%>%
  distinct()

Med$Next<-as.numeric(Med$Next)
summary(Med$Next)

rm(Long, LongSum, Med, Smoke, CatSmoke)
