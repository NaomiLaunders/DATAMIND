#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pull in the files from my master clean version and clean for this analysis
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~Libraries
library(dplyr)
library(lubridate)

#Clear environment

rm(list = ls(all.names = TRUE))

#Load full clean data

load("DatamindCPRD/Data/FinalCohort.Rdata")

####Other useful variables####
CPRD$TimeSinceDiag<-CPRD$end-CPRD$FirstSMIDate
CPRD$TimeSinceDiag<-as.numeric(CPRD$TimeSinceDiag/365.25)
summary(CPRD$TimeSinceDiag)

CPRD$AgeAtStart<-year(CPRD$start)-CPRD$yob
CPRD$AgeAtStart<-as.numeric(CPRD$AgeAtStart)
summary(CPRD$AgeAtStart)

CPRD$AgeAtEnd<-year(CPRD$end)-CPRD$yob
CPRD$AgeAtEnd<-as.numeric(CPRD$AgeAtEnd)
summary(CPRD$AgeAtEnd)


####Sort deaths and check dates####
summary(CPRD$deathdate)
#Died after exit but within follow up - didnt die in this study
length(which(CPRD$deathdate>"2018-03-31"))
CPRD$died[CPRD$deathdate>"2018-03-31"]<-0
CPRD$deathdate[CPRD$deathdate>"2018-03-31"]<-NA

#Died after 100, but within follow up - didnt die in this study
length(which(CPRD$deathdate>CPRD$Date100))
CPRD$died[CPRD$deathdate>CPRD$Date100]<-0
CPRD$deathdate[CPRD$deathdate>CPRD$Date100]<-NA

#Who died after end date (185 patients)
length(which(CPRD$deathdate>CPRD$end))
length(which(CPRD$deathdate>CPRD$end&CPRD$end==CPRD$tod))

#If died more than 6 months after tod they didnt die
length(which(CPRD$died==1&is.na(CPRD$deathdate)))
CPRD$DeathTime<-CPRD$deathdate-CPRD$end
length(which(CPRD$DeathTime>182&CPRD$deathdate>CPRD$end))
CPRD$died[CPRD$DeathTime>182&CPRD$deathdate>CPRD$end]<-0
CPRD$deathdate[CPRD$DeathTime>182&CPRD$deathdate>CPRD$end]<-NA

#If died within 6 months, end date become date of death
length(which(CPRD$DeathTime<=182&CPRD$DeathTime>0))
died<-subset(CPRD,CPRD$DeathTime<=182&CPRD$DeathTime>0)
died<-select(died, patid, deathdate, end, tod, lcd, Date100, start, crd, Date18, source, DeathTime)
length(which(died$DeathTime>30))
CPRD$end[CPRD$DeathTime<=182&CPRD$DeathTime>0&!is.na(CPRD$DeathTime)]<-CPRD$deathdate[CPRD$DeathTime<=182&CPRD$DeathTime>0&!is.na(CPRD$DeathTime)]

save(CPRD, file="DatamindCPRD/Data/CleanFinalCohort.Rdata")
