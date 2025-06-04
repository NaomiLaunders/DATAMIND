#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pull in the files from my master clean version and clean for this analysis
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~Libraries
library(tidyverse)

#Clear environment

rm(list = ls(all.names = TRUE))

#Load full clean data

load("FullCleanData/AllSMIComorbs.Rdata")

####Check dates####

#Aurum and Gold date fields are named differently so change them to be in the same field

AllSMI$deathdate<-if_else(!is.na(AllSMI$deathdate), AllSMI$deathdate, AllSMI$cprd_ddate)
AllSMI$tod<-if_else(!is.na(AllSMI$tod), AllSMI$tod, AllSMI$regenddate)
AllSMI$crd<-if_else(!is.na(AllSMI$regstartdate), AllSMI$regstartdate, AllSMI$crd)
AllSMI<-select(AllSMI, -regstartdate, -regenddate, -cprd_ddate, -regend)

#deathdate: Check all patients that died have a date of death. They do. 
length(which(AllSMI$died==1&is.na(AllSMI$deathdate)))

####Generate new fields/missing data####

#dob: Create a date of birth from year. Set at 01/01 so technically including patients who were diagnosed in the year that they turned 18.
AllSMI$dob<-paste0("01/01/", AllSMI$yob)
AllSMI$dob<-as.Date(AllSMI$dob, "%d/%m/%Y")

#Date turned 18 and 100
AllSMI$Date18<-AllSMI$dob+years(18)
AllSMI$Date100<-AllSMI$dob+years(100)

CPRD<-AllSMI

CPRD<-select(CPRD, patid, pracid, tod, crd, lcd, gender, yob, dob, Date18, Date100, deathdate, source, died, region, ethn_white, ethn_missing, 
             Broad, Ch_Total, Elix_Total, Lau_Total, Alcohol, Drugs, Date_Alcohol, Date_Drugs) 


#Age at death
CPRD$AgeAtDeath<-CPRD$deathdate-CPRD$dob
CPRD$AgeAtDeath<-as.numeric(CPRD$AgeAtDeath)
CPRD$AgeAtDeath<-CPRD$AgeAtDeath/365.25

CPRD$NoComorb<-0
CPRD$NoComorb[CPRD$Lau_Total==0]<-1
CPRD$OneComorb<-0
CPRD$OneComorb[CPRD$Lau_Total==1]<-1
CPRD$AtLeastOneComorb<-0
CPRD$AtLeastOneComorb[CPRD$Lau_Total>=1]<-1
CPRD$MoreComorb<-0
CPRD$MoreComorb[CPRD$Lau_Total>1]<-1

CPRD$NoComorb<-as.factor(CPRD$NoComorb)
CPRD$OneComorb<-as.factor(CPRD$OneComorb)
CPRD$MoreComorb<-as.factor(CPRD$MoreComorb)
CPRD$AtLeastOneComorb<-as.factor(CPRD$AtLeastOneComorb)

#Ethnicity
table(CPRD$ethn_missing, CPRD$ethn_white, useNA="ifany")
CPRD$ethn_missing[is.na(CPRD$ethn_missing)]<-"Missing"

CPRD$ethn_missing<-as.factor(CPRD$ethn_missing)
CPRD$ethn_white<-as.factor(CPRD$ethn_white)

#Save clean file
save(CPRD, file = "DatamindCPRD/Data/CleanCPRD.Rdata")

#Clear environment

rm(list = ls(all.names = TRUE))
