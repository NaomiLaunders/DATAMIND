rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(tableone)
library(stringr)

####Select only observations that are BP - GOLD####
rm(list = ls(all.names = TRUE))
BPGold<-read.table("CodeLists/Blood pressure/NaomiBPGoldFinalScreenAndValue.txt", header=TRUE, colClasses=c(medcode="character"))
BPGold<-select(BPGold, medcode, QOF, term=readterm, type)

#Count how many observation files you have in the Gold data file which start with the word observation and end with .txt
Observation_files <- list.files(path="/Raw data/Observations/NewGold/Simp", pattern = ("Obs.+\\.txt$"))
Observation_files<-paste0("", Observation_files)

#Create a table to store your results
PatBPGold<-data.frame()

#Select those that might have values - dont include record 20 as that is referral.

for (i in (1:19)) {
  FileName<-Observation_files[i]
  load(FileName)
  PatObs<-select(PatObs, patid, eventdate, sysdate, medcode, enttype, adid)
  PatObs<-merge(x=PatObs, y=BPGold, by="medcode", all.x=FALSE, all.y=FALSE)
  PatBPGold<-rbind(PatBPGold, PatObs)
  print(Observation_files[i])
}
rm(PatObs)

#If eventdate is <=01/01/1900 take sysdate
PatBPGold$eventdate<-as.Date(PatBPGold$eventdate, format="%d/%m/%Y")
PatBPGold$sysdate<-as.Date(PatBPGold$sysdate, format="%d/%m/%Y")

length(which(is.na(PatBPGold$eventdate)))
length(which(is.na(PatBPGold$sysdate)))

length(which(PatBPGold$eventdate<="1900/01/01"))
length(which(PatBPGold$sysdate<="1900/01/01"))

PatBPGold$eventdate[PatBPGold$eventdate<="1900/01/01"]<-NA
PatBPGold$sysdate[PatBPGold$sysdate<="1900/01/01"]<-NA

PatBPGold$eventdate[is.na(PatBPGold$eventdate)]<-PatBPGold$sysdate[is.na(PatBPGold$eventdate)]


save(PatBPGold, file="VariableExtracts/CPRD2018/MergedObs/GoldBPvalue.Rdata")

#Now look at test files
PatBPGold1<-data.frame()

for (i in (21:38)) {
  FileName<-Observation_files[i]
  load(FileName)
  PatObs<-select(PatObs, patid, eventdate, sysdate, medcode, enttype, data1, data2, data3, data4, data5, data6, data7, data8)
  PatObs<-merge(x=PatObs, y=BPGold, by="medcode", all.x=FALSE, all.y=FALSE)
  PatBPGold1<-rbind(PatBPGold1, PatObs)
  print(Observation_files[i])
}
rm(PatObs)

PatBPGold1$eventdate<-as.Date(PatBPGold1$eventdate, format="%d/%m/%Y")
PatBPGold1$sysdate<-as.Date(PatBPGold1$sysdate, format="%d/%m/%Y")

length(which(is.na(PatBPGold1$eventdate)))
length(which(is.na(PatBPGold1$sysdate)))

length(which(PatBPGold1$eventdate<="1900/01/01"))
length(which(PatBPGold1$sysdate<="1900/01/01"))

PatBPGold1$eventdate[PatBPGold1$eventdate<="1900/01/01"]<-NA
PatBPGold1$sysdate[PatBPGold1$sysdate<="1900/01/01"]<-NA

PatBPGold1$eventdate[is.na(PatBPGold1$eventdate)]<-PatBPGold1$sysdate[is.na(PatBPGold1$eventdate)]

save(PatBPGold1, file="VariableExtracts/CPRD2018/MergedObs/GoldBPtest.Rdata")

#####LEAVE THE VALUE FILES FOR WHEN WE WANT TO LOOK AT THEM LATER!!!####
#Now look at all the files as these may denote screening
#select all records
PatBPGold2<-data.frame()

for (i in (length(Observation_files))) {
  FileName<-Observation_files[i]
  load(FileName)
  PatObs<-select(PatObs, patid, eventdate, sysdate, medcode)
  PatObs<-merge(x=PatObs, y=BPGold, by="medcode", all.x=FALSE, all.y=FALSE)
  PatBPGold2<-rbind(PatBPGold2, PatObs)
  print(Observation_files[i])
}
rm(PatObs)


PatBPGold2$eventdate<-as.Date(PatBPGold2$eventdate, format="%d/%m/%Y")
PatBPGold2$sysdate<-as.Date(PatBPGold2$sysdate, format="%d/%m/%Y")

length(which(is.na(PatBPGold2$eventdate)))
length(which(is.na(PatBPGold2$sysdate)))

length(which(PatBPGold2$eventdate<="1900/01/01"))
length(which(PatBPGold2$sysdate<="1900/01/01"))

PatBPGold2$eventdate[PatBPGold2$eventdate<="1900/01/01"]<-NA
PatBPGold2$sysdate[PatBPGold2$sysdate<="1900/01/01"]<-NA

PatBPGold2$eventdate[is.na(PatBPGold2$eventdate)]<-PatBPGold2$sysdate[is.na(PatBPGold2$eventdate)]

save(PatBPGold2, file="VariableExtracts/CPRD2018/MergedObs/GoldBPscreen.Rdata")

####Select only observations that are SMI - AURUM####
rm(list = ls(all.names = TRUE))

NaomiAurum<-read.table("CodeLists/Blood pressure/NaomiBPAurumFinalScreenAndValue.txt", header=TRUE, colClasses=c(MedCodeId="character"))

BPAurum<-select(NaomiAurum, medcodeid = MedCodeId, QOF, Term, type)

#Count how many observation files you have in the Gold data file which start with the word observation and end with .txt
Observation_files <- list.files(path="/Raw data/Observations/NewAurum/Simp", pattern = ("Obs.+\\.txt$"))
Observation_files<-paste0("", Observation_files)

#Create a table to store your results
PatBPGold<-data.frame()

#Loop through observations
PatBPAurum<-data.frame()

for (i in (length(Observation_files))) {
  FileName<-Observation_files[i]
  load(FileName)
  PatObs<-select(PatObs, patid, medcodeid, obsdate, enterdate, value, numunitid, obstypeid, numrangelow, numrangehigh)
  PatObs<-merge(x=PatObs, y=BPAurum, by="medcodeid", all.x=FALSE, all.y=FALSE)
  PatBPAurum<-rbind(PatBPAurum, PatObs)
  print(Observation_files[i])
}
rm(PatObs)

PatBPAurum<-rename(PatBPAurum, medcode=medcodeid, sysdate=enterdate, eventdate=obsdate)

#If eventdate is <=01/01/1900 take sysdate

PatBPAurum$eventdate<-as.Date(PatBPAurum$eventdate, format="%d/%m/%Y")
PatBPAurum$sysdate<-as.Date(PatBPAurum$sysdate, format="%d/%m/%Y")

length(which(is.na(PatBPAurum$eventdate)))
length(which(is.na(PatBPAurum$sysdate)))

length(which(PatBPAurum$eventdate<="1900/01/01"))
length(which(PatBPAurum$sysdate<="1900/01/01"))

PatBPAurum$eventdate[PatBPAurum$eventdate<="1900/01/01"]<-NA
PatBPAurum$sysdate[PatBPAurum$sysdate<="1900/01/01"]<-NA

PatBPAurum$eventdate[is.na(PatBPAurum$eventdate)]<-PatBPAurum$sysdate[is.na(PatBPAurum$eventdate)]
save(PatBPAurum, file="VariableExtracts/CPRD2018/MergedObs/AurumBPvalue.Rdata")

save(PatBPAurum, file="VariableExtracts/CPRD2018/MergedObs/AurumBPscreen.Rdata")

####Merge Aurum and Gold files and bind to full SMI cohort
load("VariableExtracts/CPRD2018/MergedObs/GoldBPscreen.Rdata")
PatBPAurum<-select(PatBPAurum, -value, -numunitid, -obstypeid, -numrangelow, -numrangehigh)
PatBPAurum<-rename(PatBPAurum, term=Term)

PatBPAll<-rbind(PatBPAurum, PatBPGold2)

save(PatBPAll, file="VariableExtracts/CPRD2018/MergedObs/BPScreen.Rdata")