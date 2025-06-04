#Double check all lists
rm(list = ls(all.names = TRUE))
library(haven)
library(psych)
library(tidyverse)
library(tableone)
library(stringr)

#Load code lists

####Dictionaries####
AurumMed<-read_dta("Dictionary Browsers/CPRD_CodeBrowser_201904_Aurum/CPRDAurumMedical.dta")
AurumMed$term<-tolower(AurumMed$term)
AurumMed$term<-gsub('"', '', AurumMed$term)

GoldMed<-read.csv("Naomi/CodeLists/medicalGold.csv", header=FALSE, fill=TRUE)
GoldMed$V7<-tolower(GoldMed$V7)
GoldMed$V1<-as.character(GoldMed$V1)

####Blood pressure####
BPGold<-read.table("CodeLists/Blood pressure/NaomiBPGoldFinalScreenAndValue.txt",  sep="\t", header=TRUE, colClasses=c(medcode="character"))
BPAurum<-read.table("CodeLists/Blood pressure/NaomiBPAurumFinalScreenAndValue.txt",  sep="\t", header=TRUE, colClasses=c(MedCodeId="character"))

BPAurum$Term<-tolower(BPAurum$Term)
BPGold$readterm<-tolower(BPGold$readterm)

CheckAurum<-subset(BPGold, !(readterm %in% BPAurum$Term) & !(readcode %in% BPAurum$CleansedReadCode))
CheckGold<-subset(BPAurum, !(Term %in% BPGold$readterm))
CheckGold<-subset(CheckGold, CleansedReadCode!="")
CheckGold<-subset(CheckGold, !(CleansedReadCode %in% BPGold$readcode))

#None to add

table(BPGold$type, useNA="ifany")
table(BPAurum$type, useNA="ifany")
BPGold$type[BPGold$type=="Screening"]<-"screening"
BPAurum$type[BPAurum$type=="Screening"]<-"screening"

GoldScreen<-subset(BPGold, type=="screening")
AurumScreen<-subset(BPAurum, type=="screening")
CheckAurum<-subset(GoldScreen, !(readterm %in% AurumScreen$Term))
CheckGold<-subset(AurumScreen, !(Term %in% GoldScreen$readterm))
CheckGold<-subset(CheckGold, CleansedReadCode!="")
CheckGold<-subset(CheckGold, !(CleansedReadCode %in% GoldScreen$readcode))

BPGold<-BPGold%>%
  mutate(type = case_when(medcode=="12948" | medcode=="21826"~ "value",
                          TRUE ~ type))

#Which are included as values
BPGold$type[grepl("invit", BPGold$readterm)]<-"value"
BPGold$type[grepl("letter", BPGold$readterm)]<-"value"
BPGold$type[grepl("admin", BPGold$readterm)]<-"value"
BPGold$type[grepl("H/O", BPGold$readterm)]<-"value"
BPGold$type[grepl("Refuse", BPGold$readterm)]<-"value"
BPGold$type[grepl("deleted", BPGold$readterm)]<-"value"
BPGold$type[grepl("stopped", BPGold$readterm)]<-"value"
BPGold$type[grepl("decline", BPGold$readterm)]<-"value"
BPGold$type[grepl("hypertension monit", BPGold$readterm)]<-"value"
BPGold$type[grepl("hypertension", BPGold$readterm) & grepl("review", BPGold$readterm)]<-"value"
BPGold$type[grepl("hypertension", BPGold$readterm) & grepl("clinic", BPGold$readterm)]<-"value"
BPGold$type[grepl("hypertension self-management", BPGold$readterm)]<-"value"

BPAurum$type[grepl("invit", BPAurum$Term)]<-"value"
BPAurum$type[grepl("letter", BPAurum$Term)]<-"value"
BPAurum$type[grepl("admin", BPAurum$Term)]<-"value"
BPAurum$type[grepl("H/O", BPAurum$Term)]<-"value"
BPAurum$type[grepl("Refuse", BPAurum$Term)]<-"value"
BPAurum$type[grepl("deleted", BPAurum$Term)]<-"value"
BPAurum$type[grepl("stopped", BPAurum$Term)]<-"value"
BPAurum$type[grepl("decline", BPAurum$Term)]<-"value"
BPAurum$type[grepl("hypertension monit", BPAurum$Term)]<-"value"
BPAurum$type[grepl("hypertension", BPAurum$Term) & grepl("review", BPAurum$Term)]<-"value"
BPAurum$type[grepl("hypertension self-management", BPAurum$Term)]<-"value"

#Also check refer, monitor, hypertension, 

GoldScreen<-subset(BPGold, type=="screening")
AurumScreen<-subset(BPAurum, type=="screening")
GoldVal<-subset(BPGold, type!="screening")
AurumVal<-subset(BPAurum, type!="screening")
BPGold$type[BPGold$readterm=="seen in hypertension clinic"]<-"screening"
BPGold$type[BPGold$readterm=="seen in hypertension clinic"]<-"screening"

TestGold<-select(BPGold, RC=readcode, GoldTerm=readterm, GoldQOF=QOF, GoldType=type)
TestAurum<-select(BPAurum, RC=CleansedReadCode, AurumTerm=Term, AurumQOF=QOF, AurumType=type)

Test<-merge(x=TestGold, y=TestAurum, by="RC", all.x=FALSE, all.y=FALSE)
Test<-subset(Test, GoldQOF!=AurumQOF | GoldType!=AurumType)

write.table(BPGold, file = "CodeLists/Blood pressure/NaomiBPGoldFinalScreenAndValue.txt",  col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")
write.table(BPAurum, file = "CodeLists/Blood pressure/NaomiBPAurumFinalScreenAndValue.txt",  col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")

####Final checks for inconsistency####
CheckGold<-BPGold%>%ungroup()%>%select(readcode, readterm, QOF, type)
CheckAurum<-BPAurum%>%ungroup()%>%select(readcode=CleansedReadCode, readterm=Term, QOF, type)
Check<-rbind(CheckGold, CheckAurum)

Check<-Check%>%
  group_by(readcode)%>%
  mutate(n=n())%>%
  group_by(readcode, QOF)%>%
  mutate(QOFCheck=n())%>%
  group_by(readcode, type)%>%
  mutate(TypeCheck=n())

Check<-subset(Check, n==2 & (QOFCheck!=2 |TypeCheck!=2))

