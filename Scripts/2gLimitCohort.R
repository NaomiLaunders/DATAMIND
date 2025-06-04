rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(dplyr)
library(lubridate)
library(tidyverse)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/CPRD_FinalFull.Rdata")

####load SMI observation file####

load("VariableExtracts/CPRD2018/MergedObs/SMIObs.Rdata")

####does everyone have a code as per my list####
CountCheck<-SMIObs%>%
  group_by(patid)%>%
  summarise(nSMIDiagAll=n())

length(which(CPRD$patid %in% SMIObs$patid))
SMIObs<-subset(SMIObs, patid %in% CPRD$patid)

####Do they all have an SMI Diagnosis date####
length(which(CPRD$patid %in% SMIObs$patid))
Check<-subset(CPRD, is.na(FirstSMIDate))

#Create a table to store the results of exclusions so can easily create a flow chart for pubs
Consort<-data.frame("Count" = length(CPRD$patid), "Description" = c("Original"))

##1/Diagnosis date
CPRD<-subset(CPRD, !is.na(FirstSMIDate))
Consort1<-data.frame("Count" =length(CPRD$patid), "Description" = c("Missing diagnosis date"))
Consort<- rbind(Consort,Consort1 )

####Do they all have a diagnosis before end of follow up?####
##2/Diagnosis after end of follow up
CPRD<-subset(CPRD, FirstSMIDate<=as.Date("2018-03-31"))
Consort1<-data.frame("Count" =length(CPRD$patid), "Description" = c("Diagnosed after the end of follow up"))
Consort<- rbind(Consort,Consort1 )

##3/Under 18 at SMI diagnosis
CPRD<-subset(CPRD, year(FirstSMIDate)-yob>=18)
Consort1<-data.frame("Count" =length(CPRD$patid), "Description" = c("Diagnosed under 18"))
Consort<- rbind(Consort,Consort1 )

##4/#Drop those over 100 at index
CPRD<-subset(CPRD,(year(FirstSMIDate)-yob<=100))
Consort1<- data.frame("Count" =length(CPRD$patid), "Description" = c("Over 100 at SMI diagnosis"))
Consort<- rbind(Consort,Consort1)

##5/ End before index date
CPRD<-subset(CPRD,(end>FirstSMIDate))
Consort1<-data.frame("Count" =length(CPRD$patid), "Description" = c("Exit before SMI diagnosis"))
Consort<- rbind(Consort,Consort1)

##5/ Less than 1 year follow up
CPRD$StudyFU<-CPRD$end-CPRD$start
CPRD$StudyFU<-as.numeric(CPRD$StudyFU/365.25)
summary(CPRD$StudyFU)
CPRD<-subset(CPRD, StudyFU>=1)
Consort1<- data.frame("Count" =length(CPRD$patid), "Description" = c("Less than 1 year's follow up"))
Consort<- rbind(Consort,Consort1)

##6/ Missing region
CPRD<-subset(CPRD, !(is.na(CPRD$region)))
Consort1<- data.frame("Count" =length(CPRD$patid), "Description" = c("Missing region"))
Consort<- rbind(Consort,Consort1)

####Save numbers for flow chart####

#Generate consort
Consort
Consort<-Consort %>%
  mutate (Lag = lag(Count), Diff = Lag-Count)
Consort<-select(Consort, -Lag)
save(Consort, file = "DatamindCPRD/Outputs/FlowChart.Rdata")
write.csv(Consort, file = "DatamindCPRD/Outputs/FlowChart.csv")

####What is the time period of the study?####
summary(CPRD$start)
summary(CPRD$end)
summary(CPRD$FirstSMIDate)
table(CPRD$source)

#Save clean file
save(CPRD, file = "DatamindCPRD/Data/FinalCohort.Rdata")

#Clear environment

rm(list = ls(all.names = TRUE))
