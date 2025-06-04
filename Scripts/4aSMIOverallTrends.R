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

####Load SMI cohort####

load("DatamindCPRD/Data/FinalCohortActiveFactors.Rdata")

#####Combined measure of screening - value, category or screening####

####Study start and end####
Dates<-select(CPRD, patid, start, end, FirstSMIDate)

####BMI Combined####
load("VariableExtracts/CPRD2018/MergedObs/BMIAll.Rdata")
BMIScreen<-subset(PatBMIAll, type=="screening")
BMIScreen<-select(BMIScreen, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/BMIValue.Rdata")
BMIValue<-select(BMIValue, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/BMICalc.Rdata")
BMICalc<-select(BMICalc, patid, eventdate)
BMI<-distinct(rbind(BMIScreen, BMIValue, BMICalc))

BMI<-subset(BMI, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BMI<-subset(BMI, patid %in% CPRD$patid)
BMI<-merge(x=BMI, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMI<-subset(BMI, eventdate>=start & eventdate<=end)

save(BMI, file = "DatamindCPRD/Data/BMIAnalysis.Rdata")

####Lipids all measures####
load("VariableExtracts/CPRD2018/MergedObs/LipidScreen.Rdata")
LipidScreen<-subset(PatLipidAll, type=="screening")
LipidScreen<-select(LipidScreen, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/FinalLipids.Rdata")
LipidValue<-select(FinalLipids, patid, eventdate)

Lipid<-distinct(rbind(LipidScreen, LipidValue))

Lipid<-subset(Lipid, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Lipid<-subset(Lipid, patid %in% CPRD$patid)
Lipid<-merge(x=Lipid, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipid<-subset(Lipid, eventdate>=start & eventdate<=end)

save(Lipid, file = "DatamindCPRD/Data/LipidAnalysis.Rdata")

####BP All####
load("VariableExtracts/CPRD2018/MergedObs/BPScreen.Rdata")
BPScreen<-subset(PatBPAll, type=="screening")
BPScreen<-select(BPScreen, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/BPValues.Rdata")
BPValue<-select(PatBPAll, patid, eventdate)

BP<-distinct(rbind(BPScreen, BPValue))
BP<-subset(BP, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BP<-subset(BP, patid %in% CPRD$patid)
BP<-merge(x=BP, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BP<-subset(BP, eventdate>=start & eventdate<=end)

save(BP, file = "DatamindCPRD/Data/BPAnalysis.Rdata")

####Glucose All####
load("VariableExtracts/CPRD2018/MergedObs/GlucoseScreen.Rdata")
GlucScreen<-subset(PatGlucoseAll, type=="screening")
GlucScreen<-subset(GlucScreen, !grepl("urine", term))
GlucScreen<-select(GlucScreen, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/GlucoseValues.Rdata")
FinalGlucose<-subset(FinalGlucose, !(is.na(Hba1c) & is.na(Glucose)))
GlucValue<-select(FinalGlucose, patid, eventdate)

Gluc<-distinct(rbind(GlucScreen, GlucValue))
Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Gluc<-subset(Gluc, patid %in% CPRD$patid)
Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)

save(Gluc, file = "DatamindCPRD/Data/GlucAnalysis.Rdata")

####Smoking####
load("VariableExtracts/CPRD2018/MergedObs/SmokeScreen.Rdata")
SmokeScreen<-subset(PatSmokeAll, type=="screening")
SmokeScreen<-select(SmokeScreen, patid, eventdate)
load("VariableExtracts/CPRD2018/MergedObs/Smoke_FinalCat.Rdata")
SmokeCat<-ungroup(CatSmoke)
SmokeCat<-select(SmokeCat, patid, eventdate)

Smoke<-distinct(rbind(SmokeScreen, SmokeCat))
Smoke<-subset(Smoke, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Smoke<-subset(Smoke, patid %in% CPRD$patid)
Smoke<-merge(x=Smoke, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Smoke<-subset(Smoke, eventdate>=start & eventdate<=end)

save(Smoke, file = "DatamindCPRD/Data/SmokeAnalysis.Rdata")

####Alcohol####
load("VariableExtracts/CPRD2018/MergedObs/AlcScreen.Rdata")

Alc<-subset(PatAlcAll, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Alc<-subset(Alc, patid %in% CPRD$patid)
Alc<-distinct(select(Alc, patid, eventdate))
Alc<-merge(x=Alc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Alc<-subset(Alc, eventdate>=start & eventdate<=end)

save(Alc, file = "DatamindCPRD/Data/AlcAnalysis.Rdata")

####Ever screened for at least one thing####
CPRD$AllBMI<-0
CPRD$AllBMI[CPRD$patid %in% BMI$patid]<-1

CPRD$AllLipid<-0
CPRD$AllLipid[CPRD$patid %in% Lipid$patid]<-1

CPRD$AllGluc<-0
CPRD$AllGluc[CPRD$patid %in% Gluc$patid]<-1

CPRD$AllBP<-0
CPRD$AllBP[CPRD$patid %in% BP$patid]<-1

CPRD$AllSmoke<-0
CPRD$AllSmoke[CPRD$patid %in% Smoke$patid]<-1

CPRD$AllAlc<-0
CPRD$AllAlc[CPRD$patid %in% Alc$patid]<-1

CPRD$AnyScreen<-0
CPRD$AnyScreen[CPRD$AllBMI==1|CPRD$AllLipid==1|CPRD$AllGluc==1|CPRD$AllBP==1|
                 CPRD$AllSmoke==1|CPRD$AllAlc==1]<-1

CPRD$AllScreen<-0
CPRD$AllScreen[CPRD$AllBMI==1&CPRD$AllLipid==1&CPRD$AllGluc==1&CPRD$AllBP==1&
                   CPRD$AllSmoke==1&CPRD$AllAlc==1]<-1


which(names(CPRD)== 'AllBMI')
which(names(CPRD)== 'AnyScreen')
which(names(CPRD)== 'AllScreen')

CPRD[,c(73, 167:173)]<-lapply(CPRD[,c(73, 167:173)], factor)

#QOF
CPRD$OtherQOF<-0
CPRD$OtherQOF[!is.na(CPRD$FirstSpecificQOF)]<-1
CPRD$OtherQOF<-as.factor(CPRD$OtherQOF)

#Exemptions
load("VariableExtracts/CPRD2018/MergedObs/SMIEx.Rdata")
PatExAll<-subset(PatExAll, patid %in% CPRD$patid & Group=="MH")
Exempt<-subset(PatExAll, eventdate>="2000-04-01" & eventdate<="2018-03-31")

CPRD$EverMHExAll<-0
CPRD$EverMHExAll[CPRD$patid %in% Exempt$patid]<-1


save(CPRD, file = "DatamindCPRD/Data/MainAnalysis.Rdata")

#####Basic analysis#####
table(CPRD$AnyScreen)
prop.table(table(CPRD$AnyScreen))*100

table(CPRD$EverMHEx)
table(CPRD$EverMHEx[CPRD$AnyScreen==0])
prop.table(table(CPRD$EverMHEx[CPRD$AnyScreen==0]))*100

table(CPRD$EverMHEx[CPRD$AnyScreen==1])
prop.table(table(CPRD$EverMHEx[CPRD$AnyScreen==1]))*100

table(CPRD$AllScreen)
prop.table(table(CPRD$AllScreen))*100

table(CPRD$EverMHEx[CPRD$AllScreen==0])
prop.table(table(CPRD$EverMHEx[CPRD$AllScreen==0]))*100

table(CPRD$EverMHEx[CPRD$AllScreen==1])
prop.table(table(CPRD$EverMHEx[CPRD$AllScreen==1]))*100

CPRD$FU<-as.numeric((CPRD$end-CPRD$start)/365.25)

MyVars<-c("AgeAtSMI", "AgeAtStart", "AgeAtEnd", "gender", "ethn_missing", "country", "LastDiag", "died", "AgeAtDeath", "EverMHExAll", "FU", "OtherQOF", "AnyPres", "AnyScreen", "AllScreen", "AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc")
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table5_Screening.csv")

####BY year####
####Set active CPRD cohort####
CPRDActive<-CPRD%>%
  select(patid, Act0001, Act0102, Act0203, Act0304, Act0405, Act0506, Act0607, Act0708, Act0809, Act0910, Act1011, Act1112, Act1213, Act1314, Act1415, Act1516, Act1617, Act1718)

CPRDActive[,c(2:19)]<-lapply(CPRDActive[,c(2:19)], as.character)
CPRDActive[,c(2:19)]<-lapply(CPRDActive[,c(2:19)], as.numeric)
 
#BMI
AnyBMI<-BMI%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnyBMI<-AnyBMI%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllBMI<-merge(x=CPRDActive, y=AnyBMI, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
BMIVar<-vars_select(names(AnyBMI), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(BMIVar))) {
  AllBMI[paste0(BMIVar[i])]<-if_else(AllBMI[paste0(ActVar[i])]==0, 0, if_else(AllBMI[paste0(ActVar[i])]==1 & AllBMI[paste0(BMIVar[i])]>0, 1, 0))
  AllBMI[paste0(BMIVar[i])][is.na(AllBMI[paste0(BMIVar[i])])]<-0
}

#Calculate totals

FinalBMI<-AllBMI%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllBMI%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

BMIFin<-cbind(AllActive, FinalBMI)

BMIFin$Anyprop<-BMIFin$AnyCount/BMIFin$ActCount*100
BMIFin$Year<-str_replace(BMIFin$ActYear, "Act", "")
BMIFin$Measure<-"BMI"

#Lipid
AnyLipid<-Lipid%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnyLipid<-AnyLipid%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllLipid<-merge(x=CPRDActive, y=AnyLipid, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
LipidVar<-vars_select(names(AnyLipid), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(LipidVar))) {
  AllLipid[paste0(LipidVar[i])]<-if_else(AllLipid[paste0(ActVar[i])]==0, 0, if_else(AllLipid[paste0(ActVar[i])]==1 & AllLipid[paste0(LipidVar[i])]>0, 1, 0))
  AllLipid[paste0(LipidVar[i])][is.na(AllLipid[paste0(LipidVar[i])])]<-0
}

#Calculate totals

FinalLipid<-AllLipid%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllLipid%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

LipidFin<-cbind(AllActive, FinalLipid)

LipidFin$Anyprop<-LipidFin$AnyCount/LipidFin$ActCount*100
LipidFin$Year<-str_replace(LipidFin$ActYear, "Act", "")
LipidFin$Measure<-"Cholesterol"
#BP
AnyBP<-BP%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnyBP<-AnyBP%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllBP<-merge(x=CPRDActive, y=AnyBP, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
BPVar<-vars_select(names(AnyBP), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(BPVar))) {
  AllBP[paste0(BPVar[i])]<-if_else(AllBP[paste0(ActVar[i])]==0, 0, if_else(AllBP[paste0(ActVar[i])]==1 & AllBP[paste0(BPVar[i])]>0, 1, 0))
  AllBP[paste0(BPVar[i])][is.na(AllBP[paste0(BPVar[i])])]<-0
 }

#Calculate totals

FinalBP<-AllBP%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllBP%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

BPFin<-cbind(AllActive, FinalBP)

BPFin$Anyprop<-BPFin$AnyCount/BPFin$ActCount*100
BPFin$Year<-str_replace(BPFin$ActYear, "Act", "")
BPFin$Measure<-"Blood pressure"
#Gluc
AnyGluc<-Gluc%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnyGluc<-AnyGluc%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllGluc<-merge(x=CPRDActive, y=AnyGluc, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
GlucVar<-vars_select(names(AnyGluc), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(GlucVar))) {
  AllGluc[paste0(GlucVar[i])]<-if_else(AllGluc[paste0(ActVar[i])]==0, 0, if_else(AllGluc[paste0(ActVar[i])]==1 & AllGluc[paste0(GlucVar[i])]>0, 1, 0))
  AllGluc[paste0(GlucVar[i])][is.na(AllGluc[paste0(GlucVar[i])])]<-0
}

#Calculate totals

FinalGluc<-AllGluc%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllGluc%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

GlucFin<-cbind(AllActive, FinalGluc)

GlucFin$Anyprop<-GlucFin$AnyCount/GlucFin$ActCount*100
GlucFin$Year<-str_replace(GlucFin$ActYear, "Act", "")
GlucFin$Measure<-"Glucose"
#Smoke
AnySmoke<-Smoke%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnySmoke<-AnySmoke%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllSmoke<-merge(x=CPRDActive, y=AnySmoke, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
SmokeVar<-vars_select(names(AnySmoke), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(SmokeVar))) {
  AllSmoke[paste0(SmokeVar[i])]<-if_else(AllSmoke[paste0(ActVar[i])]==0, 0, if_else(AllSmoke[paste0(ActVar[i])]==1 & AllSmoke[paste0(SmokeVar[i])]>0, 1, 0))
  AllSmoke[paste0(SmokeVar[i])][is.na(AllSmoke[paste0(SmokeVar[i])])]<-0
}

#Calculate totals

FinalSmoke<-AllSmoke%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllSmoke%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

SmokeFin<-cbind(AllActive, FinalSmoke)

SmokeFin$Anyprop<-SmokeFin$AnyCount/SmokeFin$ActCount*100
SmokeFin$Year<-str_replace(SmokeFin$ActYear, "Act", "")
SmokeFin$Measure<-"Smoking"

#Alc
AnyAlc<-Alc%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))

AnyAlc<-AnyAlc%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

AllAlc<-merge(x=CPRDActive, y=AnyAlc, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
AlcVar<-vars_select(names(AnyAlc), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(AlcVar))) {
  AllAlc[paste0(AlcVar[i])]<-if_else(AllAlc[paste0(ActVar[i])]==0, 0, if_else(AllAlc[paste0(ActVar[i])]==1 & AllAlc[paste0(AlcVar[i])]>0, 1, 0))
  AllAlc[paste0(AlcVar[i])][is.na(AllAlc[paste0(AlcVar[i])])]<-0
}

#Calculate totals

FinalAlc<-AllAlc%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AllAlc%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

AlcFin<-cbind(AllActive, FinalAlc)

AlcFin$Anyprop<-AlcFin$AnyCount/AlcFin$ActCount*100
AlcFin$Year<-str_replace(AlcFin$ActYear, "Act", "")
AlcFin$Measure<-"Alcohol"

####Bind together all results####
AllResults<-rbind(BMIFin, LipidFin, BPFin, GlucFin, SmokeFin, AlcFin)

AllCI<-BinomCI(AllResults$AnyCount, AllResults$ActCount)*100
AllResults<-cbind(AllResults, AllCI)

cbPalette <- c("#999999", "#E69F00","#009E73", "#F0E442", "#0072B2",  "#CC79A7")

AllResults<-select(AllResults, Measure, Year, Count=AnyCount, Denom=ActCount, Prevalence=est, lower=lwr.ci, upper=upr.ci)
AllResults$Year<-case_when(AllResults$Year=="0001" ~ "2000-2001",
                           AllResults$Year=="0102" ~ "2001-2002",
                           AllResults$Year=="0203" ~ "2002-2003",
                           AllResults$Year=="0304" ~ "2003-2004",
                           AllResults$Year=="0405" ~ "2004-2005",
                           AllResults$Year=="0506" ~ "2005-2006",
                           AllResults$Year=="0607" ~ "2006-2007",
                           AllResults$Year=="0708" ~ "2007-2008",
                           AllResults$Year=="0809" ~ "2008-2009",
                           AllResults$Year=="0910" ~ "2009-2010",
                           AllResults$Year=="1011" ~ "2010-2011",
                           AllResults$Year=="1112" ~ "2011-2012",
                           AllResults$Year=="1213" ~ "2012-2013",
                           AllResults$Year=="1314" ~ "2013-2014",
                           AllResults$Year=="1415" ~ "2014-2015",
                           AllResults$Year=="1516" ~ "2015-2016",
                           AllResults$Year=="1617" ~ "2016-2017",
                           AllResults$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllResults, stat = "identity", aes(x=Year, y=Prevalence, color=Measure, group=Measure))+
  geom_ribbon(data = AllResults, aes(x=Year, ymin = lower, ymax = upper, group=Measure, fill=Measure), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  theme(legend.position = c(0.2, 0.85), legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving screening", fill="Cardiovascular risk factor")+
  geom_vline(xintercept=4.5, linetype='dotted')+
  annotate("text", x = 4.5, y = 10, label = "Physical health check", vjust = -0.5, angle = 90, , size=5)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 10, label = "All measures", vjust = -0.5, angle = 90, , size=5)+
  geom_vline(xintercept=14.5, linetype='dotted')+
  annotate("text", x = 14.5, y = 12, label = "Alcohol, BP and smoking", vjust = -0.5, angle = 90, size=5)
  
  
  
save(AllResults, file="DatamindCPRD/Outputs/ScreeningTrends.Rdata")
write.csv(AllResults, file="DatamindCPRD/Outputs/ScreeningTrends.csv")

ggsave("DatamindCPRD/Outputs/ScreeningTrends.pdf", width = 10,   height = 10)

####Ever received stratified#####
#Country
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="country", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningCountry.csv")

#gender
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="gender", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningGender.csv")

#Ethnicity
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="ethn_missing", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningEthnicity.csv")

#Diagnosis
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="LastDiag", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningDiagnosis.csv")

#Exempt (All)
table(CPRD$EverMHEx, CPRD$EverMHExAll)

MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="EverMHExAll", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningExempt.csv")

#OtherQOF
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="OtherQOF", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningQOF.csv")

#Age
CPRD$Age40<-0
CPRD$Age40[CPRD$AgeAtStart>=40]<-1
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="Age40", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningAge.csv")

#APs
MyVars<-c("AllBMI", "AllLipid", "AllGluc", "AllBP", "AllSmoke", "AllAlc", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  strata="AnyPres", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_ScreeningAP.csv")
