rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(lubridate)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/FinalCohortActiveFactors.Rdata")

####Active in 2017####

Active<-subset(CPRD, start<=as.Date("2017-12-31") & end>=as.Date("2017-01-01"))
table(Active$country)

####Proportion of remission codes who get a subsequent SMI code####
Remission<-subset(CPRD, SMIRemission==1)
load("VariableExtracts/CPRD2018/MergedObs/SMIObs.Rdata")
SMIObs<-subset(SMIObs, patid %in% CPRD$patid)
SMIObs<-subset(SMIObs, SMIDate<="2018-03-31")

Rem<-subset(SMIObs,  QOF=="Remission" & patid %in% CPRD$patid)
length(unique(Rem$patid))
remitpat<-subset(SMIObs, patid %in% Rem$patid)
length(unique(remitpat$patid))

#Set first remission date and limit to observations after that date
remitpat<-remitpat%>%
  mutate(RemissionDate=case_when(QOF=="Remission" ~ SMIDate,
                                 TRUE ~ as.Date(NA)))%>%
  group_by(patid)%>%
  mutate(FirstRem=min(RemissionDate, na.rm = TRUE))%>%
  subset(SMIDate>FirstRem)%>%
  subset(QOF!="Remission")

#How many are SMI codes
length(unique(remitpat$patid))

####QOF codes for SMI####
NoQOF<-subset(CPRD, SMIQOF==0)
table(NoQOF$country)
11389/13138
table(year(NoQOF$start))
table(year(CPRD$start))
table(NoQOF$pracid)
NoQOFreg<-as.data.frame(table(NoQOF$region))
NoQOFreg<-rename(NoQOFreg, region=Var1, NoQOF=Freq)
CPRDreg<-as.data.frame(table(CPRD$region))
CPRDreg<-rename(CPRDreg, region=Var1, pop=Freq)
NoQOF1<-cbind(NoQOFreg, CPRDreg)
NoQOF1$Perc<-NoQOF1$NoQOF/NoQOF1$pop*100

table(NoQOF$LastDiag)
table(NoQOF$EverAP)
7029/13138

table(NoQOF$EverBP)
2145/13138

table(NoQOF$source)

MyVars<-c("FirstDiag", "LastDiag")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="EnglandRegion", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table6_SMIvRegion.csv")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="country", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table7_SMIvCountry.csv")

#####Are those receiving screening the same as those with a value####
table(CPRD$EverBPValCat)
length(which(CPRD$EverBPScreen==0 & CPRD$EverBPValCat==1))
274/189557*100
table(CPRD$EverBPScreen)
length(which(CPRD$EverBPScreen==1 & CPRD$EverBPValCat==0))
195/189478*100

table(CPRD$EverLipidValCat)
length(which(CPRD$EverLipidScreen==0 & CPRD$EverLipidValCat==1))
121/136014*100
table(CPRD$EverLipidScreen)
length(which(CPRD$EverLipidScreen==1 & CPRD$EverLipidValCat==0))
1004/136897*100

table(CPRD$EverGlucValCat)
length(which(CPRD$EverGlucScreen==0 & CPRD$EverGlucValCat==1))
580/157522*100
table(CPRD$EverGlucScreen)
length(which(CPRD$EverGlucScreen==1 & CPRD$EverGlucValCat==0))
2515/159457*100

table(CPRD$EverSMokingCat)
length(which(CPRD$EverSmokeScreen==0 & CPRD$EverSMokingCat==1))
50/190933*100
table(CPRD$EverSmokeScreen)
length(which(CPRD$EverSmokeScreen==1 & CPRD$EverSMokingCat==0))
34/190917*100

table(CPRD$EverBMIValCat)
table(CPRD$EverBMIScreen)
120868/216136*100
length(which(CPRD$EverBMIScreen==0 & CPRD$EverBMIValCat==1))
47036/167641*100

length(which(CPRD$EverBMIScreen==0 & CPRD$EverBMIValCat==1 & CPRD$source=="Gold"))
length(which(CPRD$EverBMIScreen==1 & CPRD$source=="Gold"))
length(which(CPRD$source=="Gold"))
14577/76898*100
length(which(CPRD$EverBMIValCat==1 & CPRD$source=="Gold"))
59026/76898*100

CPRD$startyr<-year(CPRD$start)
CPRD$endyr<-year(CPRD$end)

MyVars<-c("source","startyr", "endyr", "AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "died", "AgeAtDeath", "NoComorb", "OneComorb", "MoreComorb")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="country", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table9_BasicByCOuntry.csv")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="LastDiag", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table9_BasicByDiag.csv")
