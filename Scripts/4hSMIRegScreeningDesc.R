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

Alc<-subset(Alc, patid %in% CPRD$patid)
BMI<-subset(BMI, patid %in% CPRD$patid)
Smoke<-subset(Smoke, patid %in% CPRD$patid)
Lipid<-subset(Lipid, patid %in% CPRD$patid)
Gluc<-subset(Gluc, patid %in% CPRD$patid)
BP<-subset(BP, patid %in% CPRD$patid)

####Ever had all on same day####
Smoke$cat<-"Smoke"
Alc$cat<-"Alcohol"
BMI$cat<-"BMI"
Lipid$cat<-"Lipid"
Gluc$cat<-"Glucose"
BP$cat<-"BP"

AllScreening<-rbind(Smoke, Alc, BMI, Lipid, Gluc, BP)

EverAll<-AllScreening%>%
  subset(patid %in% Smoke$patid & patid %in% BMI$patid & patid %in% Alc$patid
         & patid %in% Lipid$patid & patid %in% Gluc$patid & patid %in% BP$patid)

length(unique(EverAll$patid))  

EverAll$patdate<-paste0(EverAll$patid, "_", EverAll$eventdate)  
Smoke$patdate<-paste0(Smoke$patid, "_", Smoke$eventdate)
Alc$patdate<-paste0(Alc$patid, "_", Alc$eventdate)
BMI$patdate<-paste0(BMI$patid, "_", BMI$eventdate)
Lipid$patdate<-paste0(Lipid$patid, "_", Lipid$eventdate)
Gluc$patdate<-paste0(Gluc$patid, "_", Gluc$eventdate)
BP$patdate<-paste0(BP$patid, "_", BP$eventdate)

AllDay<-subset(EverAll, patdate %in% Smoke$patdate & patdate %in% BMI$patdate & patdate %in% Alc$patdate
               & patdate %in% Lipid$patdate & patdate %in% Gluc$patdate & patdate %in% BP$patdate)
length(unique(AllDay$patid))

#Within 30 days

#Interval of start and end of follow up for each patient
Active <- CPRD %>%
  group_by(patid)%>%
  mutate(active = iv(start, end), .keep = "unused")

#Start and end of follow up overall
bounds <- range(Active$active)
lower <- iv_start(bounds[[1]])
upper <- iv_end(bounds[[2]]) - 1L

#Every day in follow up
Time <- tibble(SevDay = seq(lower, upper, by = 1))

#Interval of every day to 7 days after - a rolling week
Time<-Time%>%
  mutate(weekend=SevDay+30,
         Week=iv(SevDay, weekend), .keep = "unused")%>%
  mutate(rowid_to_column(Time, "Row"))

#Order by patid and date of screening and create a one day event interval
Screen<-AllScreening%>%
  group_by(patid)%>%
  arrange(patid, eventdate)%>%
  mutate(Needle=iv(eventdate, eventdate+1), .keep = "unused")

#Create a uniqueid so can merge the intervals back to the main table
Screen<-rowid_to_column(Screen, var="Row")

rm(AllScreening, CPRD)

#look up where the event interval (needle) is in the rolling week interval. This creates a look up of row numbers where the intervals overlap.
locations <- iv_locate_overlaps(Screen$Needle, Time$Week)

#Merge locations back to the rolling weeks interval data, only keeping those that have an overlap 
FinalTime<-merge(x=locations, y=Time, by.x="haystack", by.y="Row", all.x=TRUE, all.y=FALSE)
rm(locations, Time)
rm(Active, Alc, allDay, BMI, bounds, BP, EverAll, Gluc, Lipid, Smoke, lower, upper)
rm(AllDay)

save(FinalTime, file="DatamindCPRD/Outputs/FinalTime.Rdata")
save(Screen, file="DatamindCPRD/Outputs/Screen.Rdata")

FinalTime<-select(FinalTime, -haystack, -weekend, -SevDay)
Screen<-select(Screen, Row, patid, cat)
#Merge back to the eventdates, only including those that have an overlap
FinalTime<-merge(x=FinalTime, y=Screen, by.x="needles", by.y="Row", all.x=TRUE, all.y=FALSE)
rm(Screen)

save(FinalTime, file="DatamindCPRD/Outputs/FinalTime.Rdata")

#Count how many of each category and summarise so have unique week interval and category for each patient
#Then group the summarised data by week to count how unique categories fall into each week 
FinalTime2<-FinalTime %>%
  group_by(patid, cat, Week)%>%
  summarise(count=n())%>%
  group_by(patid, Week)%>%
  summarise(count=n())

#Limit the intervals to those that have all six categories. Note, someone getting all six measures would be counted in multiple rolling windows, but
#I am interested if it ever happened so the number of times it happened isnt important. Counting the unique patients is all I need.
FinalTime2<-subset(FinalTime2, count==6)

#Count unique patients
length(unique(FinalTime2$patid))
