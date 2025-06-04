rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(tableone)
library(stringr)
library(lubridate)

####Load SMI cohort####

load("DatamindCPRD/Data/CleanCPRD.Rdata")

####Load BP lists####
####Aurum####
load("VariableExtracts/CPRD2018/MergedObs/AurumBPvalue.Rdata")
numunits <- read.delim("CPRD Data look up files/CPRD Aurum/NumUnit.txt")
obstype <- read.delim("CPRD Data look up files/CPRD Aurum/ObsType.txt")

BPAurumValue<-merge(x=PatBPAurum, y=numunits, by="numunitid", all.x=TRUE, all.y=FALSE)
BPAurumValue<-merge(x=BPAurumValue, y=obstype, by="obstypeid", all.x=TRUE, all.y=FALSE)
BPAurumValue<-rename(BPAurumValue, unit=Description.x, obstype=Description.y)

#Only select those with values
BPAurumValue<-subset(BPAurumValue, !is.na(value)& value!=0)
summary(BPAurumValue$value)
BPAurumValue$patdate<-paste0(BPAurumValue$patid, BPAurumValue$eventdate)

#If negative make positive

BPAurumValue$value[BPAurumValue$value<0]<-BPAurumValue$value[BPAurumValue$value<0]*-1

#Check obs types
Obs<-BPAurumValue%>%
  group_by(obstype)%>%
  summarise(Count=n())

Unit<-BPAurumValue%>%
  group_by(unit)%>%
  summarise(Count=n())

#check term if unit is systolic or diastolic
table(BPAurumValue$Term[BPAurumValue$unit=="Systolic"])
BPAurumValue$Term[BPAurumValue$unit=="Systolic" & BPAurumValue$Term=="o/e - blood pressure reading"]<- "systolic blood pressure"
table(BPAurumValue$Term[BPAurumValue$unit=="Diastolic"])
BPAurumValue$Term[BPAurumValue$unit=="Diastolic" & BPAurumValue$Term=="o/e - blood pressure reading"]<- "diastolic blood pressure"

Term<-BPAurumValue%>%
  group_by(Term)%>%
  summarise(Count=n())

TermUnit<-BPAurumValue%>%
  group_by(Term, unit)%>%
  summarise(Count=n())

#All units now make sense - most appear to be mmHg

BPAurumValue<-BPAurumValue%>%
  mutate(NewUnit=case_when(unit=="/ mmHg" | unit=="/mmHg" | unit=="mm Hg" |unit=="mm hh"|unit=="mm/Hg" |unit=="mm[Hg]"| unit=="mmHg" ~ "mmhg",
                           TRUE ~ "unknown"))

table(BPAurumValue$NewUnit)
table(BPAurumValue$NewUnit, BPAurumValue$unit)

####Label up terms####

BPAurumValue <- BPAurumValue %>%
  mutate(Systolic= case_when(grepl("systol", Term) ~ 1, 
                             TRUE ~ 0),
         Diastolic = case_when(grepl("diastol", Term) ~ 1,
                                 TRUE ~ 0),
         Mean =    case_when(Systolic==0 & Diastolic==0 ~ 1, 
                                    TRUE ~ 0))

#Which type are they from
table(BPAurumValue$type)
      
#Is the mean, really mean?
table(BPAurumValue$Term[BPAurumValue$Systolic==1])
table(BPAurumValue$Term[BPAurumValue$Diastolic==1])
table(BPAurumValue$Term[BPAurumValue$Mean==1])

BPAurumValue<-distinct(BPAurumValue)

#If multiple just take those that are sys/dia
AurumSys<-subset(BPAurumValue, Systolic==1)
AurumDia<-subset(BPAurumValue, Diastolic==1)
AurumMean<-subset(BPAurumValue, Mean==1)

AurumMean<-subset(AurumMean, !(patdate %in% AurumSys$patdate) & !(patdate %in% AurumDia$patdate))

AurumMean<-subset(AurumMean, value>20 & value<300)
AurumMean<-AurumMean%>%
  group_by(patdate)%>%
  mutate(SameDay=n(), maxval=max(value), minval=min(value))
AurumMeanKeep<-subset(AurumMean, SameDay==1)
AurumMeanKeep<-select(AurumMeanKeep, -minval, -maxval, -SameDay)

AurumMeanFix<-subset(AurumMean, SameDay>1)

AurumMean3<-AurumMeanFix%>%
  subset(SameDay==3)%>%
  group_by(patdate)%>%
  mutate(value=ave(value))%>%
  filter(row_number()==1)%>%
  select(-minval, -maxval, -SameDay)
  
AurumMeanFix<-AurumMeanFix%>%
  subset(SameDay==2)%>%
  group_by(patdate)%>%
  mutate(Mean= 0,
         Systolic = case_when(value==maxval ~ 1,
                              TRUE ~ 0),
         Diastolic = case_when(value==minval ~ 1,
                               TRUE ~ 0))%>%
  select(-minval, -maxval, -SameDay)

BPAurumValue<-rbind(AurumSys, AurumDia, AurumMeanKeep, AurumMean3, AurumMeanFix)

#Change to long format
BPAurumLong<-BPAurumValue%>%
  pivot_longer(17:19, names_to = "Label", values_to = "Count")%>%
  subset(Count>0)%>%
  select(-Count)

table(BPAurumLong$Label, useNA="ifany")

tapply(BPAurumLong$value, BPAurumLong$Label, summary)

rm(Term, BPAurumValue, PatBPAurum, Unit, obstype, numunits, Obs)
rm(AurumSys, AurumDia, AurumMeanKeep, AurumMean3, AurumMeanFix)
rm(AurumMean, AurumMulti, Check, TermUnit)

####GOLD#####
load("VariableExtracts/CPRD2018/MergedObs/GoldBPvalue.Rdata")
load("VariableExtracts/CPRD2018/MergedObs/GoldBPtest.Rdata")
entity <- read.delim("CPRD Data look up files/CPRD GOLD/entity.txt")
opr <- read.delim("CPRD Data look up files/CPRD GOLD/TXTFILES/OPR.txt")
sum <- read.delim("CPRD Data look up files/CPRD GOLD/TXTFILES/SUM.txt")
tqu <- read.delim("CPRD Data look up files/CPRD GOLD/TXTFILES/TQU.txt")

Ent<-select(entity, enttype, description)

Ad1<-read.table("Naomi/Raw data/Observations/GoldAdditional/18_288R_GOLD_Extract_Additional_001.txt", header=TRUE, fill=TRUE, sep="\t", colClasses = c("patid"="character"))
Ad2<-read.table("Naomi/Raw data/Observations/GoldAdditional/18_288R_GOLD_Extract_Additional_002.txt", header=TRUE, fill=TRUE, sep="\t", colClasses = c("patid"="character"))
Ad3<-read.table("Naomi/Raw data/Observations/GoldAdditional/18_288R_GOLD_Extract_Additional_003.txt", header=TRUE, fill=TRUE, sep="\t", colClasses = c("patid"="character"))

Additional<-rbind(Ad1, Ad2, Ad3)
rm(Ad1, Ad2, Ad3)

#How many of those with additional info are for blood pressure
Additional<-subset(Additional, patid %in% CPRD$patid)
Additional<-merge(x=Additional, y=Ent, by="enttype", all.x=TRUE, all.y=FALSE)
table(Additional$description)
BPAdd<-subset(Additional, description=="Blood pressure")
length(which(BPAdd$patid %in% PatBPGold$patid))

BPAdd$patad<-paste0(BPAdd$patid, "-", BPAdd$adid)
PatBPGold$patad<-paste0(PatBPGold$patid, "-", PatBPGold$adid)
length(which(BPAdd$patad %in% PatBPGold$patad))

#WHat is the code for those we havent got
Check<-subset(BPAdd, !(patad %in% PatBPGold$patad))

PatBPCheck<-data.frame()

for (i in (length(Observation_files))) {
  FileName<-Observation_files[i]
  load(FileName)
  PatObs$patad<-paste0(PatObs$patid, "-", PatObs$adid)
  PatObs<-subset(PatObs, PatObs$patad %in% Check$patad)
  PatBPCheck<-rbind(PatBPCheck, PatObs)
  print(Observation_files[i])
}
rm(PatObs)

##Now check the medcodes
GoldMed<-read.csv("Naomi/CodeLists/medicalGold.csv", header=FALSE, fill=TRUE)
GoldMed$V7<-tolower(GoldMed$V7)
GoldMed$V1<-as.character(GoldMed$V1)
GoldMed<-rename(GoldMed, medcode=V1)

PatBPCheck<-merge(x=PatBPCheck, y=GoldMed, by="medcode", all.x=TRUE, all.y=FALSE)

#Decide to keep all of them
Date<-select(PatBPGold, patad, eventdate)
Date2<-select(PatBPCheck, patad, eventdate)
Date2$eventdate<-as.Date(Date2$eventdate, format="%d/%m/%Y")
length(which(is.na(Date2$eventdate)))
Date<-rbind(Date, Date2)

#Add in date
BPAdd<-merge(x=BPAdd, y = Date, by="patad", all.x=TRUE, all.y=FALSE)
BPAdd<-merge(x=BPAdd, y=PatBPGold, by="patad", all.x=TRUE, all.y=FALSE)
BPAdd<-select(BPAdd, medcode, patid=patid.x, eventdate=eventdate.x, sysdate, Diastolic = data1, Systolic = data2, QOF, term, enttype=enttype.x, description)
BPAdd<-subset(BPAdd, !((is.na(Systolic)|Systolic<=0) & (is.na(Diastolic)|Diastolic<=0)))

#Try the test table
BPGoldValue<-merge(x=PatBPGold1, y=Ent, by="enttype", all.x=TRUE, all.y=FALSE)
EntCheck<-subset(entity, enttype %in% BPGoldValue$enttype)
EntCheckAdd<-subset(entity, enttype %in% BPAdd$enttype)

#Look at categories
TQUVals<-subset(BPGoldValue, enttype==392)
TQUVals<-subset(TQUVals, data1>0)
TQUVals<-merge(x=TQUVals, y=tqu, by.x="data1", by.y="Code", all.x=TRUE, all.y=FALSE)

TQUVals2<-subset(BPGoldValue, enttype==411|enttype==288 )
TQUVals2<-subset(TQUVals2, data4>0)
TQUVals2<-merge(x=TQUVals2, y=tqu, by.x="data4", by.y="Code", all.x=TRUE, all.y=FALSE)
TQUVals<-rbind(TQUVals, TQUVals2)

save(TQUVals, file="VariableExtracts/CPRD2018/MergedObs/GoldBPCat.Rdata")

#Now look at values
GoldValue<-subset(BPGoldValue, enttype==411|enttype==288)
GoldValue<-merge(x=GoldValue, y=opr, by.x="data1", by.y="Code", all.x=TRUE, all.y=FALSE)
table(GoldValue$Operator)
GoldValue<-merge(x=GoldValue, y=sum, by.x="data3", by.y="Code", all.x=TRUE, all.y=FALSE)
table(GoldValue$Specimen.Unit.of.Measure)

GoldValue<-rename(GoldValue, unit=Specimen.Unit.of.Measure, value=data2)
GoldValue<-subset(GoldValue, !is.na(value)& value!=0)

#All are postural drop and therefore not blood pressure reading. Remove.

#Back to additional table

#Change to long format
BPAddLong<-BPAdd%>%
  pivot_longer(5:6, names_to = "Label", values_to = "value")%>%
  subset(value>0)

table(BPAddLong$Label, BPAddLong$description)

BPAddLong$value<-as.numeric(BPAddLong$value)
tapply(BPAddLong$value, BPAddLong$Label, summary)

rm(BPGoldValue, opr, sum, tqu, TQUVals, TQUVals2, GoldValue, Ent, EntCheck, Check, entity, PatBPGold1)
rm(Additional, BPAdd, Date, Date2, EntCheckAdd, GoldMed, PatBPCheck, PatBPGold)
rm(TermCheck)

####Merge Gold and Aurum and clean####
BPAurumLong<-select(BPAurumLong, medcode, patid, eventdate, value, QOF, Term, unit, patdate, NewUnit, Label)
BPAddLong$unit<-NA
BPAddLong$NewUnit<-NA
BPAddLong$patdate<-paste0(BPAddLong$patid, BPAddLong$eventdate)
GoldValueLong<-select(BPAddLong, medcode, patid, eventdate, value, QOF, Term=term, unit, patdate, NewUnit, Label)

BP<-rbind(GoldValueLong, BPAurumLong)
BP<-distinct(BP)
BP$patdatelabel<-paste0(BP$patid, BP$eventdate, BP$Label)
length(unique(BP$patdatelabel))

#Drop medcode  and QOF as likely the only difference for some
BP<-select(BP, -medcode, -QOF, -Term)
BP<-distinct(BP)
length(unique(BP$patdatelabel))

save(BP, file = "VariableExtracts/CPRD2018/MergedObs/GoldandAurumBP.Rdata")

####Find those with multiple records of the same type on the same day####
Multi<-BP%>%
  group_by(patdatelabel)%>%
  mutate(SameDay=n())%>%
  subset(SameDay>1)
length(unique(Multi$patdatelabel))

Single<-BP%>%
  group_by(patdatelabel)%>%
  mutate(SameDay=n())%>%
  subset(SameDay==1)
length(unique(Single$patdatelabel))

Multi<-Multi%>%
  group_by(patdatelabel)%>%
  mutate(maxval=max(value), minval=min(value), ave=mean(value), diff=maxval-minval)%>%
  ungroup%>%
  mutate(InRange = case_when(Label=="Systolic" & value<=269.5 & value>=52.2 ~ 1,
                             Label=="Diastolic" & value<=167.2 & value>=27.9 ~ 1,
                             Label=="Mean" & value<=269.5 & value>=27.9 ~ 1,
                             TRUE ~ 0))%>%
  group_by(patdatelabel)%>%
  mutate(EverInRange = sum(InRange))%>%
  subset(InRange==1 | (InRange==0&EverInRange==0))%>%
  mutate(maxval=max(value), minval=min(value), ave=mean(value), diff=maxval-minval)

length(which(Multi$diff==0))
length(which(Multi$diff<=10))

NewSingle<-Multi %>%
  subset(diff<10)%>%
  group_by(patdatelabel) %>%
  mutate(value=ave) %>%
  filter(row_number()==1)%>%
  select(-maxval, -minval, -ave, -diff)
length(unique(NewSingle$patdatelabel))

Single<-rbind(Single, NewSingle)

LargeDiff<-subset(Multi, diff>10)

NewSingle<-LargeDiff %>%
  group_by(patdatelabel) %>%
  mutate(value=ave) %>%
  filter(row_number()==1)%>%
  select(-maxval, -minval, -ave, -diff)
length(unique(NewSingle$patdatelabel))

Single<-rbind(Single, NewSingle)

length(unique(Single$patdatelabel))

NewBP<-ungroup(Single)
NewBP<-select(NewBP, -InRange, -EverInRange, -NewUnit, -SameDay, -patdatelabel, -unit)

table(NewBP$Label)

FinalBP<-NewBP%>%
  group_by(patdate)%>%
  pivot_wider(names_from = Label, values_from = c(value))
length(unique(FinalBP$patdate))

#Rename
FinalBP$eventyear<-year(FinalBP$eventdate)

#Limit to valid values
table(FinalBP$Mean)
length(which(is.na(FinalBP$Systolic) & is.na(FinalBP$Diastolic)))

FinalBP$Systolic[FinalBP$Systolic<52.2|FinalBP$Systolic>269.5]<-NA
FinalBP$Diastolic[FinalBP$Diastolic<27.9|FinalBP$Diastolic>167.2]<-NA
FinalBP$Mean[FinalBP$Mean<27.9|FinalBP$Mean>269.5]<-NA

FinalBP$Drop<-0
FinalBP$Drop[is.na(FinalBP$Diastolic) & is.na(FinalBP$Systolic) & is.na(FinalBP$Mean)]<-1
table(FinalBP$Drop)

FinalBP<-FinalBP%>%
  ungroup()

FinalBP<-subset(FinalBP, Drop==0)
FinalBP<-select(FinalBP, -Drop)
save(FinalBP, file="VariableExtracts/CPRD2018/MergedObs/BPValues.Rdata")

####Find categories####
####Add in categories####
load("VariableExtracts/CPRD2018/MergedObs/GoldBPCat.Rdata")
table(TQUVals$term)   
TQUVals<-rename(TQUVals, Term = term)
TQUVals<-subset(TQUVals, Term!="o/e - bp reading:postural drop")

TQUVals$Mean<-1

#Change to long format
TQUValsLong<-TQUVals%>%
  pivot_longer(19, names_to = "Label", values_to = "Count")%>%
  subset(Count>0)%>%
  select(-Count)

TQUValsLong$patdate<-paste0(TQUValsLong$patid,TQUValsLong$eventdate)
TQUValsLong$patdatelabel<-paste0(TQUValsLong$patid,TQUValsLong$eventdate, TQUValsLong$Label)

table(TQUValsLong$Test.Qualifier)
TQUValsLong$Test.Qualifier[TQUValsLong$Test.Qualifier=="Abnormal"]<-"Abnormal"
TQUValsLong$Test.Qualifier[TQUValsLong$Test.Qualifier=="Potential Abnormal"]<-"Abnormal"
TQUValsLong$Test.Qualifier[TQUValsLong$Test.Qualifier=="Normal"]<-"Normal"

table(TQUValsLong$Term[TQUValsLong$Label=="Other"],TQUValsLong$Test.Qualifier[TQUValsLong$Label=="Other"])

TQ<-select(TQUValsLong, patid, patdate, patdatelabel, eventdate, Label, cat = Test.Qualifier)

####Gold categories part 2####

NaomiGold<-read.table("CodeLists/Blood pressure/NaomiBPGoldFinalScreenAndValue.txt", header=TRUE, fill=TRUE, sep="\t", dec = ".", colClasses=c(medcode="character"))

load("VariableExtracts/CPRD2018/MergedObs/GoldBPscreen.Rdata")
CatGold<-rename(PatBPGold2, Term=term)
CatGold<-select(CatGold, patid, medcode, Term, QOF, type, eventdate)

NaomiGold<-rename(NaomiGold, Term=readterm)

table(CatGold$Term)

CatGold <- CatGold %>%
  mutate(Systolic= case_when(grepl("systol", Term) ~ 1, 
                             TRUE ~ 0),
         Diastolic = case_when(grepl("diastol", Term) ~ 1,
                               TRUE ~ 0),
         Mean =    case_when(Systolic==0 & Diastolic==0 ~ 1, 
                             TRUE ~ 0))

NaomiCatGold<-NaomiGold%>%
  mutate(Qual = case_when(grepl("h/o", Term) | grepl("history", Term) ~ "Drop",
                          grepl("hyperten", Term) & !(grepl("screen", Term)) ~ "Diagnosis hypertension",
                          grepl("hypoten", Term) ~ "Diagnosis hypotension",
                                                      (grepl("borderline", Term)) ~ "Borderline",
                          
                          (Term=="highest brachial systolic pressure" | grepl("clinic", Term) | grepl("monitor", Term) |grepl("review", Term) | grepl("screen", Term)) ~ "Drop",
                          
                          (grepl("normal", Term) & !grepl("abnormal", Term)) | grepl("hypertension resolved", Term) ~ "Normal",
                          
                          
                          
                          (grepl("raised", Term) | grepl("high blood", Term) |
                             grepl("high", Term))~ "High",
                          
                          (grepl("low blood", Term) | grepl("reading low", Term) | grepl("very low", Term)) ~ "Low",
                          
                          (grepl("abnormal", Term)) ~ "Abnormal",
                          
                          TRUE ~ "Drop"))

CheckDrop<-subset(NaomiCatGold, Qual=="Drop") 
CheckDiag<-subset(NaomiCatGold, Qual=="Diagnosis hypertension") 
CheckDiag1<-subset(NaomiCatGold, Qual=="Diagnosis hypotension") 
table(NaomiCatGold$Term[NaomiCatGold$Qual=="High"])
table(NaomiCatGold$Term[NaomiCatGold$Qual=="Abnormal"]) 
table(NaomiCatGold$Term[NaomiCatGold$Qual=="Low"]) 
table(NaomiCatGold$Term[NaomiCatGold$Qual=="Normal"]) 
table(NaomiCatGold$Term[NaomiCatGold$Qual=="Borderline"])

NaomiCatGold<-select(NaomiCatGold, medcode, Qual)
CatGold<-merge(x=CatGold, y=NaomiCatGold, by="medcode", all.x=TRUE, all.y=FALSE)
CatGold<-subset(CatGold, Qual!="Drop")

CatGold<-CatGold%>%
  pivot_longer(7:9, names_to = "Label", values_to = "Count")%>%
  subset(Count>0)%>%
  select(-Count)

CatGold$patdate<-paste0(CatGold$patid, CatGold$eventdate)
CatGold$patdatelabel<-paste0(CatGold$patid,CatGold$eventdate, CatGold$Label)
CatGold<-select(CatGold, patid, patdate, patdatelabel, eventdate, Label, cat = Qual)
length(unique(CatGold$patdatelabel))

CatGold<-distinct(CatGold)

####Merge all Gold####
CatGold<-rbind(CatGold, TQ)

CatGold<-CatGold%>%
  group_by(patdatelabel)%>%
  mutate(SameDay = n())

CatGoldSingle<-subset(CatGold, SameDay==1)
CatGoldMulti<-subset(CatGold, SameDay>1)
table(CatGoldMulti$cat)

CatGoldMulti<-CatGoldMulti%>%
  group_by(patdatelabel)%>%
  mutate(Order = case_when(cat=="High" ~ 7,
                           cat=="Low" ~ 6,
                           cat=="Borderline" ~ 5,
                           cat == "Abnormal" ~ 4,
                           cat == "Normal" ~ 3,
                           cat == "Diagnosis - Hypertension" ~ 2,
                           cat == "Diagnosis - Hypotension" ~ 1,
                           
                           TRUE ~ 0))
CatGoldMulti<-CatGoldMulti%>%
  group_by(patdatelabel)%>%
  arrange(desc(Order)) %>%
  filter(row_number()==1)
length(unique(CatGoldMulti$patdatelabel))

CatGold<-rbind(CatGoldSingle, CatGoldMulti)

CatGold<-CatGold%>%
  group_by(patdate)%>%
  select(-Order, -SameDay, -patdatelabel)%>%
  pivot_wider(names_from = Label, names_prefix = "Cat_", values_from = c(cat))

#Need to look at Aurum categories and codes that are diagnoses
NaomiAurum<-read.table("CodeLists/Blood pressure/NaomiBPAurumFinalScreenAndValue.txt", header=TRUE, fill=TRUE, sep="\t", dec = ".", colClasses=c(MedCodeId="character"))
table(NaomiAurum$Term)

load("VariableExtracts/CPRD2018/MergedObs/AurumBPscreen.Rdata")

CatAurum <- PatBPAurum %>%
  mutate(Systolic= case_when(grepl("systol", Term) ~ 1, 
                             TRUE ~ 0),
         Diastolic = case_when(grepl("diastol", Term) ~ 1,
                               TRUE ~ 0),
         Mean =    case_when(Systolic==0 & Diastolic==0 ~ 1, 
                             TRUE ~ 0))

table(NaomiAurum$Term)
NaomiCatAurum<-NaomiAurum%>%
  mutate(Qual = 
    case_when(grepl("h/o", Term) | grepl("history", Term) ~ "Drop",
              grepl("hyperten", Term) & !(grepl("screen", Term)) ~ "Diagnosis hypertension",
              grepl("hypoten", Term) ~ "Diagnosis hypotension",
              (grepl("borderline", Term)) ~ "Borderline",
                          
                          (Term=="hypertension medication review" | Term=="highest ankle systolic pressure" | Term=="highest brachial systolic pressure" | grepl("clinic", Term) | grepl("monitor", Term) |grepl("review", Term) | grepl("screen", Term)) ~ "Drop",
                          
                          (grepl("normal", Term) & !grepl("abnormal", Term)) ~ "Normal",
              (grepl("raised", Term)| grepl("high blood", Term) |
                             grepl("high", Term))~ "High",
                          
                           grepl("low blood", Term) | grepl("reading low", Term) | grepl("very low", Term) ~ "Low",
                          
                          (grepl("abn", Term)) ~ "Abnormal",
                          
                          TRUE ~ "Drop"))

NaomiCatAurum$Qual[NaomiCatAurum$Term=="af screen using bp monitor with af detector abnormal"]<- "Abnormal"
NaomiCatAurum$Qual[NaomiCatAurum$Term=="[rfc] high/low bp"]<- "Abnormal"

CheckDrop<-subset(NaomiCatAurum, Qual=="Drop") 
CheckDiag<-subset(NaomiCatAurum, Qual=="Diagnosis hypertension") 
CheckDiag1<-subset(NaomiCatAurum, Qual=="Diagnosis hypotension") 
table(NaomiCatAurum$Term[NaomiCatAurum$Qual=="High"])
table(NaomiCatAurum$Term[NaomiCatAurum$Qual=="Abnormal"]) 
table(NaomiCatAurum$Term[NaomiCatAurum$Qual=="Low"]) 
table(NaomiCatAurum$Term[NaomiCatAurum$Qual=="Normal"]) 
table(NaomiCatAurum$Term[NaomiCatAurum$Qual=="Borderline"])

NaomiCatAurum<-select(NaomiCatAurum, medcode=MedCodeId, Qual)
CatAurum<-merge(x=CatAurum, y=NaomiCatAurum, by="medcode", all.x=TRUE, all.y=FALSE)
CatAurum<-subset(CatAurum, Qual!="Drop")

CatAurum<-CatAurum%>%
  pivot_longer(13:15, names_to = "Label", values_to = "Count")%>%
  subset(Count>0)%>%
  select(-Count)

CatAurum$patdate<-paste0(CatAurum$patid, CatAurum$eventdate)
CatAurum$patdatelabel<-paste0(CatAurum$patid,CatAurum$eventdate, CatAurum$Label)
CatAurum<-select(CatAurum, patid, patdate, patdatelabel, eventdate, Label, cat = Qual)
length(unique(CatAurum$patdatelabel))

CatAurum<-distinct(CatAurum)

CatAurum<-CatAurum%>%
  group_by(patdatelabel)%>%
  mutate(SameDay = n())

CatAurumSingle<-subset(CatAurum, SameDay==1)
CatAurumMulti<-subset(CatAurum, SameDay>1)
table(CatAurumMulti$cat)

CatAurumMulti<-CatAurumMulti%>%
  group_by(patdatelabel)%>%
  mutate(Order = case_when(cat=="High" ~ 7,
                           cat =="Low" ~ 6,
                           cat=="Borderline" ~ 5,
                           cat == "Abnormal" ~ 4,
                           cat == "Normal" ~ 3,
                           cat == "Diagnosis - Hypertension" ~ 2,
                           cat == "Diagnosis - Hypotension" ~ 1,
                           TRUE ~ 0))
CatAurumMulti<-CatAurumMulti%>%
  group_by(patdatelabel)%>%
  arrange(desc(Order)) %>%
  filter(row_number()==1)
length(unique(CatAurumMulti$patdatelabel))

CatAurum<-rbind(CatAurumSingle, CatAurumMulti)

CatAurum<-CatAurum%>%
  group_by(patdate)%>%
  select(-Order, -SameDay, -patdatelabel)%>%
  pivot_wider(names_from = Label, names_prefix = "Cat_", values_from = c(cat))

CatBP<-rbind(CatAurum, CatGold)
length(unique(CatBP$patdate))

AllBP<-merge(x=FinalBP, y=CatBP, by="patdate", all.x=TRUE, all.y=TRUE)

length(unique(AllBP$patdate))

AllBP<-AllBP%>%
  mutate(patid = coalesce(patid.x, patid.y),
         eventdate = coalesce(eventdate.x, eventdate.y),
         value = case_when(is.na(patid.x) ~ 0,
                           TRUE ~ 1),
         cat = case_when(is.na(patid.y) ~ 0,
                         TRUE ~ 1))%>%
  select(-patid.x, -patid.y, -eventdate.x, -eventdate.y)

save(AllBP, file = "VariableExtracts/CPRD2018/MergedObs/AllBP_FinalinclCat.Rdata")



