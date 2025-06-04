rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(lubridate)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/FinalCohortActivewSMI.Rdata")

CPRD$FirstDiag<-factor(CPRD$FirstDiag, levels=c("schizophrenia", "bipolar", "other psychosis"))
CPRD$LastDiag<-factor(CPRD$LastDiag, levels=c("schizophrenia", "bipolar", "other psychosis"))

CPRD$diagyear<-year(CPRD$FirstSMIDate)
CPRD$QOFyear<-year(CPRD$FirstQOF)

CPRD$SMIQOF<-0
CPRD$SMIQOF[!is.na(CPRD$FirstQOF)]<-1
CPRD$SMIQOF<-as.factor(CPRD$SMIQOF)

CPRD$AnyQOF<-0
CPRD$AnyQOF[!is.na(CPRD$FirstAnyQOF)]<-1
CPRD$AnyQOF<-as.factor(CPRD$AnyQOF)

#Check classes and fix
sapply(CPRD[,c(1:155)], class)

CPRD$died<-as.factor(CPRD$died)
CPRD$region<-as.factor(CPRD$region)
CPRD$Broad<-as.factor(CPRD$Broad)
CPRD$Alcohol<-as.factor(CPRD$Alcohol)
CPRD$Drugs<-as.factor(CPRD$Drugs)
CPRD$EverAP<-as.factor(CPRD$EverAP)
CPRD$EverBP<-as.factor(CPRD$EverBP)
CPRD$EverEvidenceAP<-as.factor(CPRD$EverEvidenceAP)
CPRD$EverEvidenceBP<-as.factor(CPRD$EverEvidenceBP)
CPRD$country<-as.factor(CPRD$country)
CPRD$EnglandRegion<-as.factor(CPRD$EnglandRegion)

which(names(CPRD)=="PracticeIMD")
which(names(CPRD)=="EverOtherEx")
CPRD[,30:44]<-lapply(CPRD[,30:44], factor)

which(names(CPRD)=="EverBMIScreen")
which(names(CPRD)=="Dementia")
CPRD[,48:91]<-lapply(CPRD[,48:91], factor)

which(names(CPRD)=="FirstDiagisQOF")
which(names(CPRD)=="SMIRemission")
CPRD[,149:154]<-lapply(CPRD[,149:154], factor)

which(names(CPRD)== 'Act0001')
which(names(CPRD)== 'Act1718')
CPRD[,125:142]<-lapply(CPRD[,125:142], factor)

CPRD$AnyPres<-0
CPRD$AnyPres[CPRD$EverAP==1 | CPRD$EverBP==1]<-1
CPRD$AnyPres<-as.factor(CPRD$AnyPres)

CPRD$NoPresButAPEv<-0
CPRD$NoPresButAPEv[CPRD$EverAP==0 & CPRD$EverEvidenceAP==1]<-1
CPRD$NoPresButAPEv<-as.factor(CPRD$NoPresButAPEv)

CPRD$NoPresButBPEv<-0
CPRD$NoPresButBPEv[CPRD$EverBP==0 & CPRD$EverEvidenceBP==1]<-1
CPRD$NoPresButBPEv<-as.factor(CPRD$NoPresButBPEv)

CPRD$EverBMIValCat<-0
CPRD$EverBMIValCat[CPRD$EverBMIValue==1| CPRD$EverBMICat==1]<-1
CPRD$EverBMIValCat<-as.factor(CPRD$EverBMIValCat)


CPRD$EverLipidValCat<-0
CPRD$EverLipidValCat[CPRD$EverLipidValue==1| CPRD$EverLipidCat==1]<-1
CPRD$EverLipidValCat<-as.factor(CPRD$EverLipidValCat)

CPRD$EverBPValCat<-0
CPRD$EverBPValCat[CPRD$EverBPValue==1| CPRD$EverBPCat==1]<-1
CPRD$EverBPValCat<-as.factor(CPRD$EverBPValCat)

CPRD$EverGlucValCat<-0
CPRD$EverGlucValCat[CPRD$EverGlucoseValue==1| CPRD$EverGlucoseCat==1]<-1
CPRD$EverGlucValCat<-as.factor(CPRD$EverGlucValCat)

save(CPRD, file = "DatamindCPRD/Data/FinalCohortActiveFactors.Rdata")

#Look at follow up time
ggplot(CPRD, aes(x=StudyFU)) + geom_histogram()


#Table 1 - Cohort basics
MyVars<-c("source","AgeAtSMI", "AgeAtStart", "AgeAtEnd", "StudyFU", "gender", "ethn_missing", "country", "EnglandRegion", "FullIMD", "died", "AgeAtDeath", "NoComorb", "OneComorb", "MoreComorb")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_Basic.csv")

table(CPRD$FullIMD, CPRD$PatientIMD, useNA="ifany")
length(which(!is.na(CPRD$FullIMD)))
length(which(!is.na(CPRD$FullIMD)&is.na(CPRD$PatientIMD) ))

#Table 2 - Who is active
which(names(CPRD)== 'Act0001')
which(names(CPRD)== 'Act1718')

MyVars<-names(CPRD[,125:142])
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table2_Active.csv")

#Table 3 - SMI diagnosis and antipsychotics
MyVars<-c("SMIQOF", "SMIRemission", "FirstDiag", "LastDiag", "EverSchizophrenia", "EverBipolar", "EverOther", "diagyear", "TimeSinceDiag", "QOFyear", "AnyPres", "EverAP", "EverBP", "NoPresButAPEv", "NoPresButBPEv")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table3_SMI.csv")

#Table 3 - SMI diagnosis and antipsychotics by country
MyVars<-c("SMIQOF", "SMIRemission", "FirstDiag", "LastDiag", "EverSchizophrenia", "EverBipolar", "EverOther", "diagyear", "TimeSinceDiag", "QOFyear", "AnyPres", "EverAP", "EverBP", "NoPresButAPEv", "NoPresButBPEv")

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="country", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table3a_SMI.csv")


#Table 4 - Exemptions (Note - MH exemptions are QOF only)
which(names(CPRD)== 'EverNonQOFExempt')
which(names(CPRD)== 'EverOtherEx')

MyVars<-names(CPRD[,c(33:44)])

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = TRUE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table4_Exempt.csv")

#Table 5 - Screening and QOF
which(names(CPRD)== 'EverBMIScreen')
which(names(CPRD)== 'Dementia')

MyVars<-names(CPRD[,c(48:91)])
MyExtra<-c("EverBMIValCat", "EverLipidValCat", "EverBPValCat", "EverGlucValCat")
MyVars<-append(MyVars, MyExtra )

Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table5_Screening.csv")

#Table 5 - Screening and QOF
which(names(CPRD)== 'EverBMIScreen')
which(names(CPRD)== 'Dementia')

MyVars<-names(CPRD[,c(48:91)])
MyExtra<-c("EverBMIValCat", "EverLipidValCat", "EverBPValCat", "EverGlucValCat")
MyVars<-append(MyVars, MyExtra )

Table1<-CreateTableOne(vars=MyVars,  strata="source", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table5_ScreeningBySource.csv")

#Table 5 - Screening and country
Table1<-CreateTableOne(vars=MyVars,  strata="country", data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, showAllLevels = TRUE)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE, showAllLevels = FALSE)
write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table5_ScreeningbyCountry.csv")
