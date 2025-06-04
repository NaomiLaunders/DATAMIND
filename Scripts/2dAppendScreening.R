rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/CPRD_SMIEx.Rdata")

####Study start and end####
CPRD$start<-pmax(CPRD$crd, CPRD$Date18, CPRD$FirstSMIDate, as.Date("2000-04-01"), na.rm=TRUE)
CPRD$end<-pmin(as.Date("2018-03-31"), CPRD$tod, CPRD$deathdate, CPRD$lcd, CPRD$Date100, na.rm=TRUE)
Dates<-select(CPRD, patid, start, end, FirstSMIDate)

####BMI screening####
load("VariableExtracts/CPRD2018/MergedObs/BMIAll.Rdata")
BMIScreen<-subset(PatBMIAll, type=="screening")
BMI<-subset(BMIScreen, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BMI<-subset(BMI, patid %in% CPRD$patid)
BMI<-merge(x=BMI, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMI<-subset(BMI, eventdate>=start & eventdate<=end)
BMIQOF<-subset(BMI, QOF==1)

CPRD$EverBMIScreen<-0
CPRD$EverBMIScreen[CPRD$patid %in% BMI$patid]<-1
CPRD$EverBMIQOF<-0
CPRD$EverBMIQOF[CPRD$patid %in% BMIQOF$patid]<-1

prop.table(table(CPRD$EverBMIScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBMIQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBMIQOF, CPRD$source),2)*100
rm(BMI, BMIQOF, BMIScreen)

####BMI value####
load("VariableExtracts/CPRD2018/MergedObs/BMIValue.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/BMICalc.Rdata")
BMIValue<-subset(BMIValue, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BMICalc<-subset(BMICalc, eventdate>="2000-04-01" & eventdate<="2018-03-31")

BMIValue<-subset(BMIValue, patid %in% CPRD$patid)
BMIValue<-merge(x=BMIValue, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMIValue<-subset(BMIValue, eventdate>=start & eventdate<=end)

BMICalc<-subset(BMICalc, patid %in% CPRD$patid)
BMICalc<-merge(x=BMICalc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMICalc<-subset(BMICalc, eventdate>=start & eventdate<=end)

BMIValue$patmonthyear<-paste0(BMIValue$patid, "-", month(BMIValue$eventdate), year(BMIValue$eventdate))
BMICalc$patmonthyear<-paste0(BMICalc$patid, "-", month(BMIValue$eventdate), year(BMICalc$eventdate))

#Only include calculated if value not present that month and year
BMICalc<-subset(BMICalc, !(patmonthyear %in% BMIValue$patmonthyear))
BMICalc<-rename(BMICalc, value=BMICalc)
#Merge all together
BMIAnyValue<-rbind(BMIValue, BMICalc)

CPRD$EverBMIValue<-0
CPRD$EverBMIValue[CPRD$patid %in% BMIAnyValue$patid]<-1

prop.table(table(CPRD$EverBMIValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBMIValue, CPRD$source),2)*100
rm(BMICalc, BMIValue, BMIAnyValue)

####BMI category####
load("VariableExtracts/CPRD2018/MergedObs/BMIAll.Rdata")
BMICat<-subset(PatBMIAll, !is.na(BMICat))
BMICat<-subset(BMICat, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BMICat<-subset(BMICat, patid %in% CPRD$patid)
BMICat<-merge(x=BMICat, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BMICat<-subset(BMICat, eventdate>=start & eventdate<=end)

CPRD$EverBMICat<-0
CPRD$EverBMICat[CPRD$patid %in% BMICat$patid]<-1

prop.table(table(CPRD$EverBMICat, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBMICat, CPRD$source),2)*100

rm(BMICat)

####Lipids screening####
load("VariableExtracts/CPRD2018/MergedObs/LipidScreen.Rdata")
table(PatLipidAll$type)
Lipid<-subset(PatLipidAll, type=="screening")
Lipid<-subset(Lipid, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Lipid<-subset(Lipid, patid %in% CPRD$patid)
Lipid<-merge(x=Lipid, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Lipid<-subset(Lipid, eventdate>=start & eventdate<=end)

LipidQOF<-subset(Lipid, QOF>0)

CPRD$EverLipidScreen<-0
CPRD$EverLipidScreen[CPRD$patid %in% Lipid$patid]<-1
CPRD$EverLipidQOF<-0
CPRD$EverLipidQOF[CPRD$patid %in% LipidQOF$patid]<-1

prop.table(table(CPRD$EverLipidScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverLipidQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverLipidQOF, CPRD$source),2)*100

rm(Lipid, LipidQOF, PatLipidAll)

####Lipids value####
load("VariableExtracts/CPRD2018/MergedObs/FinalLipids.Rdata")
LipidValue<-subset(FinalLipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")
LipidValue<-subset(LipidValue, patid %in% CPRD$patid)
LipidValue<-merge(x=LipidValue, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
LipidValue<-subset(LipidValue, eventdate>=start & eventdate<=end)

CPRD$EverLipidValue<-0
CPRD$EverLipidValue[CPRD$patid %in% LipidValue$patid]<-1

Tchol<-subset(LipidValue, !is.na(LipidValue$tchol))

CPRD$EverTcholValue<-0
CPRD$EverTcholValue[CPRD$patid %in% Tchol$patid]<-1

Tcholhdl<-subset(LipidValue, !is.na(LipidValue$tchdlratio))

CPRD$EverTchdlratioValue<-0
CPRD$EverTchdlratioValue[CPRD$patid %in% Tcholhdl$patid]<-1

prop.table(table(CPRD$EverLipidValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverLipidValue, CPRD$source),2)*100

prop.table(table(CPRD$EverTcholValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverTcholValue, CPRD$source),2)*100

prop.table(table(CPRD$EverTchdlratioValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverTchdlratioValue, CPRD$source),2)*100

rm(LipidValue, FinalLipids)

####Lipids category####
load("VariableExtracts/CPRD2018/MergedObs/AllLipids_FinalinclCat.Rdata")
LipidCat<-subset(AllLipids, eventdate>="2000-04-01" & eventdate<="2018-03-31")
LipidCat<-subset(LipidCat, patid %in% CPRD$patid)
LipidCat<-merge(x=LipidCat, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
LipidCat<-subset(LipidCat, eventdate>=start & eventdate<=end)

LipidCat<-subset(LipidCat, cat==1)

CPRD$EverLipidCat<-0
CPRD$EverLipidCat[CPRD$patid %in% LipidCat$patid]<-1

prop.table(table(CPRD$EverLipidCat, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverLipidCat, CPRD$source),2)*100

rm(AllLipids, LipidCat)

####BP screening####
load("VariableExtracts/CPRD2018/MergedObs/BPScreen.Rdata")
table(PatBPAll$type)
BP<-subset(PatBPAll, type=="screening")
BP<-subset(BP, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BP<-subset(BP, patid %in% CPRD$patid)
BP<-merge(x=BP, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BP<-subset(BP, eventdate>=start & eventdate<=end)

BPQOF<-subset(BP, QOF==1)

CPRD$EverBPScreen<-0
CPRD$EverBPScreen[CPRD$patid %in% BP$patid]<-1
CPRD$EverBPQOF<-0
CPRD$EverBPQOF[CPRD$patid %in% BPQOF$patid]<-1

prop.table(table(CPRD$EverBPScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBPQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBPQOF, CPRD$source),2)*100

rm(BP, BPQOF, PatBPAll)

####BP values####
load("VariableExtracts/CPRD2018/MergedObs/BPValues.Rdata")
BPValue<-subset(FinalBP, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BPValue<-subset(BPValue, patid %in% CPRD$patid)
BPValue<-merge(x=BPValue, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BPValue<-subset(BPValue, eventdate>=start & eventdate<=end)

CPRD$EverBPValue<-0
CPRD$EverBPValue[CPRD$patid %in% BPValue$patid]<-1

prop.table(table(CPRD$EverBPValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBPValue, CPRD$source),2)*100

rm(FinalBP, BPValue)

####BP category####
load("VariableExtracts/CPRD2018/MergedObs/AllBP_FinalinclCat.Rdata")

BPCat<-subset(AllBP, eventdate>="2000-04-01" & eventdate<="2018-03-31")
BPCat<-subset(BPCat, patid %in% CPRD$patid)
BPCat<-merge(x=BPCat, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
BPCat<-subset(BPCat, eventdate>=start & eventdate<=end)

BPCat<-subset(BPCat, cat==1)

CPRD$EverBPCat<-0
CPRD$EverBPCat[CPRD$patid %in% BPCat$patid]<-1

prop.table(table(CPRD$EverBPCat, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBPCat, CPRD$source),2)*100

rm(AllBP, BPCat)

####Glucose screening####
load("VariableExtracts/CPRD2018/MergedObs/GlucoseScreen.Rdata")
table(PatGlucoseAll$type)
Gluc<-subset(PatGlucoseAll, type=="screening")
Gluc<-subset(Gluc, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Gluc<-subset(Gluc, patid %in% CPRD$patid)
Gluc<-merge(x=Gluc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Gluc<-subset(Gluc, eventdate>=start & eventdate<=end)

GlucQOF<-subset(Gluc, QOF==1)

CPRD$EverGlucScreen<-0
CPRD$EverGlucScreen[CPRD$patid %in% Gluc$patid]<-1
CPRD$EverGlucQOF<-0
CPRD$EverGlucQOF[CPRD$patid %in% GlucQOF$patid]<-1

prop.table(table(CPRD$EverGlucScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverGlucQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverGlucQOF, CPRD$source),2)*100

rm(Gluc, GlucQOF, PatGlucoseAll)

####Glucose values####
load("VariableExtracts/CPRD2018/MergedObs/GlucoseValues.Rdata")
GlucoseValue<-subset(FinalGlucose, eventdate>="2000-04-01" & eventdate<="2018-03-31")
GlucoseValue<-subset(GlucoseValue, patid %in% CPRD$patid)
GlucoseValue<-merge(x=GlucoseValue, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
GlucoseValue<-subset(GlucoseValue, eventdate>=start & eventdate<=end)

CPRD$EverGlucoseValue<-0
CPRD$EverGlucoseValue[CPRD$patid %in% GlucoseValue$patid]<-1

BloodGlucose<-subset(GlucoseValue, !is.na(GlucoseValue$Glucose))

CPRD$EverBloodGlucoseValue<-0
CPRD$EverBloodGlucoseValue[CPRD$patid %in% BloodGlucose$patid]<-1

Hba1c<-subset(GlucoseValue, !is.na(GlucoseValue$Hba1c))

CPRD$EverHba1cValue<-0
CPRD$EverHba1cValue[CPRD$patid %in% Hba1c$patid]<-1

prop.table(table(CPRD$EverGlucoseValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverGlucoseValue, CPRD$source),2)*100

prop.table(table(CPRD$EverBloodGlucoseValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverBloodGlucoseValue, CPRD$source),2)*100

prop.table(table(CPRD$EverHba1cValue, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverHba1cValue, CPRD$source),2)*100

rm(FinalGlucose, GlucoseValue)

####Glucose category####
load("VariableExtracts/CPRD2018/MergedObs/AllGlucose_FinalinclCat.Rdata")

GlucoseCat<-subset(AllGlucose, eventdate>="2000-04-01" & eventdate<="2018-03-31")
GlucoseCat<-subset(GlucoseCat, patid %in% CPRD$patid)
GlucoseCat<-merge(x=GlucoseCat, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
GlucoseCat<-subset(GlucoseCat, eventdate>=start & eventdate<=end)
GlucoseCat<-subset(GlucoseCat, cat==1)


CPRD$EverGlucoseCat<-0
CPRD$EverGlucoseCat[CPRD$patid %in% GlucoseCat$patid]<-1

prop.table(table(CPRD$EverGlucoseCat, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverGlucoseCat, CPRD$source),2)*100

rm(AllGlucose, GlucoseCat)

####Smoking####
load("VariableExtracts/CPRD2018/MergedObs/SmokeScreen.Rdata")
table(PatSmokeAll$type)
Smoke<-subset(PatSmokeAll, type=="screening")
Smoke<-subset(Smoke, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Smoke<-subset(Smoke, patid %in% CPRD$patid)
Smoke<-merge(x=Smoke, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Smoke<-subset(Smoke, eventdate>=start & eventdate<=end)
SmokeQOF<-subset(Smoke, QOF==1)

CPRD$EverSmokeScreen<-0
CPRD$EverSmokeScreen[CPRD$patid %in% Smoke$patid]<-1
CPRD$EverSmokeQOF<-0
CPRD$EverSmokeQOF[CPRD$patid %in% SmokeQOF$patid]<-1

prop.table(table(CPRD$EverSmokeScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverSmokeQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverSmokeQOF, CPRD$source),2)*100

rm(PatSmokeAll, Smoke, SmokeQOF)

####Smoking categories####
load("VariableExtracts/CPRD2018/MergedObs/Smoke_FinalCat.Rdata")

SmokingCat<-subset(CatSmoke, eventdate>="2000-04-01" & eventdate<="2018-03-31")
SmokingCat<-subset(SmokingCat, patid %in% CPRD$patid)
SmokingCat<-merge(x=SmokingCat, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
SmokingCat<-subset(SmokingCat, eventdate>=start & eventdate<=end)

CPRD$EverSMokingCat<-0
CPRD$EverSMokingCat[CPRD$patid %in% SmokingCat$patid]<-1

prop.table(table(CPRD$EverSMokingCat, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverSMokingCat, CPRD$source),2)*100

####Alcohol####
load("VariableExtracts/CPRD2018/MergedObs/AlcScreen.Rdata")

Alc<-subset(PatAlcAll, eventdate>="2000-04-01" & eventdate<="2018-03-31")
Alc<-subset(Alc, patid %in% CPRD$patid)
Alc<-merge(x=Alc, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
Alc<-subset(Alc, eventdate>=start & eventdate<=end)

AlcQOF<-subset(Alc, QOF==1)

CPRD$EverAlcScreen<-0
CPRD$EverAlcScreen[CPRD$patid %in% Alc$patid]<-1
CPRD$EverAlcQOF<-0
CPRD$EverAlcQOF[CPRD$patid %in% AlcQOF$patid]<-1

prop.table(table(CPRD$EverAlcScreen, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverAlcQOF, CPRD$FirstDiag),2)*100
prop.table(table(CPRD$EverAlcQOF, CPRD$source),2)*100

####Any screening####
CPRD$AnyScreen<-0
CPRD$AnyScreen[CPRD$EverSmokeScreen==1|CPRD$EverGlucScreen==1|CPRD$EverBPScreen==1|CPRD$EverBMIScreen==1|
                 CPRD$EverLipidScreen==1|CPRD$EverAlcScreen==1]<-1

CPRD$QOFScreen<-0
CPRD$QOFScreen[CPRD$EverGlucScreen==1|CPRD$EverBPScreen==1|CPRD$EverBMIScreen==1|
                 CPRD$EverLipidScreen==1|CPRD$EverAlcScreen==1]<-1

save(CPRD, file = "DatamindCPRD/Data/CPRD_Screen.Rdata")
