rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(lubridate)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/FinalCohortActiveFactors.Rdata")

which(names(CPRD)== 'EverBMIScreen')
which(names(CPRD)== 'EverSMokingCat')

MyVars<-names(CPRD[,c(48:70)])
MyExtra<-c("EverBMIValCat", "EverLipidValCat", "EverBPValCat", "EverGlucValCat")
MyVars<-append(MyVars, MyExtra )

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE, strata="country")
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table7_ScreeningByCountry.csv")

Graph<-select(CPRD, EverBMIValue, EverTcholValue, EverTchdlratioValue, EverBPValue, EverBloodGlucoseValue, EverHba1cValue, EverSMokingCat, country)
Vars<-names(Graph[,c(1:7)])

Graph<-Graph%>%
  group_by(country)%>%
  mutate_at(Vars, as.character)%>%
  mutate_at(Vars, as.numeric)%>%
  summarise(across(c(1:7), sum), Total=n())
  
Graph<-Graph%>%
  pivot_longer(cols=2:8, names_to="cat", values_to="Count")

Graph$Percent<-Graph$Count/Graph$Total*100

Graph$Cat[Graph$cat=="EverBMIValue"]<-"BMI"
Graph$Cat[Graph$cat=="EverTcholValue"]<-"Total cholesterol"
Graph$Cat[Graph$cat=="EverTchdlratioValue"]<-"TC:HDL"
Graph$Cat[Graph$cat=="EverBPValue"]<-"Blood pressure"
Graph$Cat[Graph$cat=="EverBloodGlucoseValue"]<-"Blood glucose"
Graph$Cat[Graph$cat=="EverHba1cValue"]<-"HbA1c"
Graph$Cat[Graph$cat=="EverSMokingCat"]<-"Smoking"

Graph$Cat<-factor(Graph$Cat, levels=c("BMI", "Total cholesterol", "TC:HDL", "Blood pressure", "Blood glucose",
                                      "HbA1c","Smoking" ))

ggplot()+
  geom_bar(data = Graph, stat = "identity", position="dodge", aes(x=Cat, y=Percent, group=country, fill=country))+
  theme_classic()+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))+
  theme(axis.text.x = element_text(size = 14, angle=90))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Risk factors", y = "Percentage of patients with value", fill="Nation")

####Longitudinal trends by country####
####BMI####
load("VariableExtracts/CPRD2018/MergedObs/BMIValue.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/BMICalc.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/BMICat.Rdata")

BMICalc<-select(BMICalc, patid, eventdate)
BMIValue<-select(BMIValue, patid, eventdate)
BMICat<-select(BMICat, patid, eventdate)

BMI<-rbind(BMICalc, BMIValue, BMICat)

BMI<-subset(BMI, patid %in% CPRD$patid)
length(unique(BMI$patid))

BMI<-subset(BMI, eventdate>="2000-04-01" & eventdate<="2018-03-31")
length(unique(BMI$patid))

Dates<-select(CPRD, patid, start, end)
BMI<-merge(x=BMI, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMI<-subset(BMI, eventdate>=start & eventdate<=end)
length(unique(BMI$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-BMI%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="BMI")

Longitudinal<-LongSum

####Total chol####

load("VariableExtracts/CPRD2018/MergedObs/AllLipids_FinalinclCat.Rdata")

AllLipids<-subset(AllLipids, patid %in% CPRD$patid)

length(unique(AllLipids$patid))
Lipids<-subset(AllLipids, !is.na(tchol) | !is.na(Cat_totalcholesterol))
length(unique(Lipids$patid))

Lipids<-subset(Lipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
Lipids<-merge(x=Lipids, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipids<-subset(Lipids, eventdate>=start & eventdate<=end)
length(unique(Lipids$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Lipids%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="Total cholesterol")

Longitudinal<-rbind(Longitudinal, LongSum)

####TC:HDL####
Lipids<-subset(AllLipids, !is.na(tchdlratio)| !is.na(Cat_tchdlratio))
length(unique(Lipids$patid))

Lipids<-subset(Lipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
Lipids<-merge(x=Lipids, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipids<-subset(Lipids, eventdate>=start & eventdate<=end)
length(unique(Lipids$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Lipids%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="TC:HDL")

Longitudinal<-rbind(Longitudinal, LongSum)

####Blood pressure####

load("VariableExtracts/CPRD2018/MergedObs/AllBP_FinalinclCat.Rdata")

AllBP<-subset(AllBP, patid %in% CPRD$patid)

length(unique(AllBP$patid))

BP<-subset(AllBP, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
BP<-merge(x=BP, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BP<-subset(BP, eventdate>=start & eventdate<=end)

#Take first record and see how many have in 6 months, 1 year and 2 years
Long<-BP%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
        upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="Blood pressure")

Longitudinal<-rbind(Longitudinal, LongSum)

####Blood glucose####

load("VariableExtracts/CPRD2018/MergedObs/AllGlucose_FinalinclCat.Rdata")

AllGlucose<-subset(AllGlucose, patid %in% CPRD$patid)

length(unique(AllGlucose$patid))

Gluc<-subset(AllGlucose, !is.na(Glucose) | !is.na(Cat_Glucose))

Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)
length(unique(Gluc$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Gluc%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="Blood glucose")

Longitudinal<-rbind(Longitudinal, LongSum)

####HbA1c####
Gluc<-subset(AllGlucose, !is.na(Hba1c) | !is.na(Cat_Hba1c))
Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)
length(unique(Gluc$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Gluc%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="HbA1c")

Longitudinal<-rbind(Longitudinal, LongSum)

####Smoking####

load("VariableExtracts/CPRD2018/MergedObs/Smoke_FinalCat.Rdata")

CatSmoke<-subset(CatSmoke, patid %in% CPRD$patid)

length(unique(CatSmoke$patid))

Smoke<-subset(CatSmoke, eventdate>="2000-04-01" & eventdate<="2018-03-31")

Dates<-select(CPRD, patid, start, end)
Smoke<-merge(x=Smoke, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Smoke<-subset(Smoke, eventdate>=start & eventdate<=end)
length(unique(Smoke$patid))

#Take first record and see how many have in 6 months, 1 year and 2 years

Long<-Smoke%>%
  group_by(patid)%>%
  mutate(First=min(eventdate),
         upto2yr=case_when(eventdate-First>7 & eventdate-First<=730.5 ~ 1, TRUE ~ 0))

#Summarise
LongSum<-Long%>%
  group_by(patid)%>%
  summarise(upto2yr=sum(upto2yr))%>%
  mutate(upto2yr=case_when(upto2yr>0~1, TRUE ~0), Cat="Smoking")

Longitudinal<-rbind(Longitudinal, LongSum)

####Final longitudinal####
Country<-select(CPRD, patid, country)
FinalLong<-merge(x=Longitudinal, y=Country, by="patid", all.x=FALSE, all.y=FALSE)
FinalLong<-FinalLong%>%
  group_by(Cat, country)%>%
  summarise(upto2yr=sum(upto2yr))
CountryDen<-CPRD%>%
  group_by(country)%>%
  summarise(Total=n())

FinalLong<-merge(x=FinalLong, y=CountryDen, by="country", all.x=FALSE, all.y=FALSE)
FinalLong$Percent<-FinalLong$upto2yr/FinalLong$Total*100  

FinalLong$Cat<-factor(FinalLong$Cat, levels=c("BMI", "Total cholesterol", "TC:HDL", "Blood pressure", "Blood glucose",
                                      "HbA1c","Smoking" ))

ggplot()+
  geom_bar(data = FinalLong, stat = "identity", position="dodge", aes(x=Cat, y=Percent, group=country, fill=country))+
  theme_classic()+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))+
  theme(axis.text.x = element_text(size = 14, angle=90))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Risk factors", y = "Percentage of patients with 2 values within 2 years", fill="Nation")
