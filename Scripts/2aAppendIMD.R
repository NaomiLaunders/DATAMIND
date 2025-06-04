####Deal with IMD####

rm(list = ls(all.names = TRUE))
library(tidyverse)
library(lubridate)
library(haven)
load("DatamindCPRD/Data/CleanCPRD.Rdata")

#Load IMD
AurumPracIMD<-read.table("Data/Linked data/Aurum-linked/practice_imd_18_288R.txt", header=TRUE, quote="", fill=TRUE, sep="\t", colClasses = ("pracid"="character"))
GoldPracIMD<-read.table("Data/Linked data/GOLD-linked/practice_imd_18_288R.txt", header=TRUE, quote="", fill=TRUE, sep="\t", colClasses = ("pracid"="character"))
AurumPatIMD<-read.table("Data/Linked data/Aurum-linked/patient_imd2015_18_288R.txt", header=TRUE, quote="", fill=TRUE, sep="\t", colClasses = ("patid"="character"))
GoldPatIMD<-read.table("Data/Linked data/GOLD-linked/patient_imd2015_18_288R.txt", header=TRUE, quote="", fill=TRUE, sep="\t", colClasses = ("patid"="character"))

#Check IMD by practice for the patients IMD
Test<-AurumPatIMD
Test$imd2015_5[Test$imd2015_5==""]<-"missing"
Test$imd2015_5[Test$imd2015_5=="1"]<-"One"
Test$imd2015_5[Test$imd2015_5=="2"]<-"Two"
Test$imd2015_5[Test$imd2015_5=="3"]<-"Three"
Test$imd2015_5[Test$imd2015_5=="4"]<-"Four"
Test$imd2015_5[Test$imd2015_5=="5"]<-"Five"

Test$Value<-1
Test<-Test%>%
  pivot_wider(names_from="imd2015_5", values_from="Value")%>%
  group_by(pracid)%>%
  summarise(Total= n(), One=sum(One, na.rm =TRUE ), Two=sum(Two, na.rm =TRUE), Three=sum(Three, na.rm =TRUE), Four=sum(Four,na.rm =TRUE), Five=sum(Five, na.rm =TRUE), Missing=sum(missing, na.rm =TRUE))

Prop<-Test%>%
  group_by(pracid)%>%
  summarise(One=(One/Total*100), Two=(Two/Total*100), Three=(Three/Total*100), Four=(Four/Total*100), Five=(Five/Total*100), Missing=(Missing/Total*100))

#Sort Gold pracid to match
GoldPracIMD$pracid<-as.numeric(GoldPracIMD$pracid)

GoldPracIMD$pracid[GoldPracIMD$pracid<100&GoldPracIMD$pracid>9]<-paste0("0", GoldPracIMD$pracid[GoldPracIMD$pracid<100&GoldPracIMD$pracid>9])
GoldPracIMD$pracid[GoldPracIMD$pracid=="9"]<-"009"
GoldPracIMD$pracid[GoldPracIMD$pracid=="8"]<-"008"
GoldPracIMD$pracid[GoldPracIMD$pracid=="7"]<-"007"
GoldPracIMD$pracid[GoldPracIMD$pracid=="6"]<-"006"
GoldPracIMD$pracid[GoldPracIMD$pracid=="5"]<-"005"
GoldPracIMD$pracid[GoldPracIMD$pracid=="4"]<-"004"
GoldPracIMD$pracid[GoldPracIMD$pracid=="3"]<-"003"
GoldPracIMD$pracid[GoldPracIMD$pracid=="2"]<-"002"
GoldPracIMD$pracid[GoldPracIMD$pracid=="1"]<-"001"

GoldPracIMD$pracid<-as.character(GoldPracIMD$pracid)
#Aurum practice IMD
CPRD<-merge(x=CPRD, y=AurumPracIMD, by="pracid", all.x=TRUE, all.y=FALSE)
CPRD<-dplyr::select(CPRD, -country, -ni2017_imd_5, -s2016_imd_5, -w2014_imd_5)
CPRD<-rename(CPRD, AurumPracIMD=e2015_imd_5)
table(CPRD$AurumPracIMD, CPRD$source, useNA="ifany")

#Gold practice IMD
CPRD<-merge(x=CPRD, y=GoldPracIMD, by="pracid", all.x=TRUE, all.y=FALSE)
CPRD<-dplyr::select(CPRD, -country, -ni2017_imd_5, -s2016_imd_5, -w2014_imd_5)
CPRD<-rename(CPRD, GoldPracIMD=e2015_imd_5)
table(CPRD$GoldPracIMD, CPRD$source, useNA="ifany")

#Practice IMD
CPRD$PracticeIMD<-CPRD$GoldPracIMD
CPRD$PracticeIMD[is.na(CPRD$PracticeIMD)]<-CPRD$AurumPracIMD[is.na(CPRD$PracticeIMD)]
table(CPRD$PracticeIMD, useNA="ifany")
CPRD<-dplyr::select(CPRD, -AurumPracIMD, -GoldPracIMD)

#Aurum patient IMD
AurumPatIMD<-dplyr::select(AurumPatIMD, -pracid)
CPRD<-merge(x=CPRD, y=AurumPatIMD, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-rename(CPRD, AurumPatIMD=imd2015_5)
table(CPRD$AurumPatIMD, CPRD$source, useNA="ifany")
CPRD$AurumPatIMD[CPRD$AurumPatIMD==""]<-NA
table(CPRD$AurumPatIMD, CPRD$source, useNA="ifany")

#Gold patient IMD
GoldPatIMD<-dplyr::select(GoldPatIMD, -pracid)
CPRD<-merge(x=CPRD, y=GoldPatIMD, by="patid", all.x=TRUE, all.y=FALSE)
CPRD<-rename(CPRD, GoldPatIMD=imd2015_5)
table(CPRD$GoldPatIMD, CPRD$source, useNA="ifany")
CPRD$GoldPatIMD[CPRD$GoldPatIMD==""]<-NA
table(CPRD$GoldPatIMD, CPRD$source, useNA="ifany")

#Patient IMD
CPRD$PatientIMD<-CPRD$GoldPatIMD
CPRD$PatientIMD[is.na(CPRD$PatientIMD)]<-CPRD$AurumPatIMD[is.na(CPRD$PatientIMD)]
table(CPRD$PatientIMD, useNA="ifany")
CPRD<-dplyr::select(CPRD, -AurumPatIMD, -GoldPatIMD)

length(which(CPRD$PatientIMD!=CPRD$PracticeIMD))
table(CPRD$PatientIMD, CPRD$PracticeIMD, useNA="ifany")

#Check for one practice
length(which(CPRD$pracid=="20002"))
Test<-subset(CPRD, pracid=="20002")
Test<-dplyr::select(Test, patid, pracid, PatientIMD, PracticeIMD)

####NOTE: THIS IS ONLY IMD FOR ENGLAND!!!!~~~~~~#######
CPRD$FullIMD<-CPRD$PatientIMD
CPRD$FullIMD[is.na(CPRD$FullIMD)]<-CPRD$PracticeIMD[is.na(CPRD$FullIMD)]

table(CPRD$FullIMD, useNA="ifany")

save(CPRD, file="DatamindCPRD/Data/CPRD_IMD.Rdata")
