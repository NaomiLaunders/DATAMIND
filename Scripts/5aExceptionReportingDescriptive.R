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
library(sandwich)
library(lmtest)
####Load SMI cohort####

load("DatamindCPRD/Data/MainAnalysis.Rdata")

cbPalette <- c("#999999", "#CC79A7", "#E69F00","#009E73", "#F0E442", "#0072B2")

####Limit to QOF years####
CPRD$start<-pmax(CPRD$start, as.Date("2004-04-01"))

###Limit Scotland to March 2016###
CPRD$end<-case_when(CPRD$country=="Scotland" ~ pmin(CPRD$end, as.Date("2016-03-31")),
                    TRUE ~ CPRD$end)

CPRD[,125:142]<-lapply(CPRD[,125:142], as.character)
CPRD[,125:142]<-lapply(CPRD[,125:142], as.numeric)  

CPRD$Act1617<-case_when(CPRD$country=="Scotland" ~ 0,
                    TRUE ~ CPRD$Act1617)
CPRD$Act1718<-case_when(CPRD$country=="Scotland" ~ 0,
                        TRUE ~ CPRD$Act1718)

#Calculate active time since 2004 and limit to those who were active for at least one year
CPRD<-CPRD%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Act04 = sum(c_across(129:142)))

CPRD<-subset(CPRD, Act04>=1)

CPRD$FU<-as.numeric(CPRD$end-CPRD$start)/365.25
summary(CPRD$FU)

CPRD$AgeAtStart<-year(CPRD$start) - CPRD$yob

CPRDActive<-select(CPRD, patid, Act0405, Act0506, Act0607, Act0708, Act0809, Act0910, Act1011, Act1112, Act1213, Act1314, Act1415, Act1516, Act1617, Act1718)



####All exception reports####
load("VariableExtracts/CPRD2018/MergedObs/SMIEx.Rdata")

####Make sure Exemption obs are just those for my cohort####
length(which(CPRD$patid %in% PatExAll$patid))
PatExAll<-subset(PatExAll, patid %in% CPRD$patid)
PatExAll<-subset(PatExAll, !(is.na(eventdate)))
ExObs<-subset(PatExAll, eventdate>="2004-04-01" & eventdate<="2018-03-31")
ExObs$Type[ExObs$Type=="MH"]<-"Unsuitable"
ExObs$Type[ExObs$Type=="Not indicated/tolerated"]<-"Unsuitable"

####Group alcohol as a group####
Alc<-subset(ExObs, Group=="Other")
Alc<-select(Alc, term)
Alc<-distinct(Alc)
Alc<-subset(Alc, grepl("alc", term, ignore.case=TRUE))

ExObs$Group<-case_when(ExObs$Group=="Other" & grepl("alc", ExObs$term, ignore.case=TRUE) ~ "Alchol",
                       TRUE ~ ExObs$Group)

####Only include those that occur while the patient is eligible####
#When did it occur?
Ex<-ExObs%>%
  mutate(Year = case_when(eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
                          eventdate>="2005-04-01" & eventdate<="2006-03-31" ~ "0506",
                          eventdate>="2006-04-01" & eventdate<="2007-03-31" ~ "0607",
                          eventdate>="2007-04-01" & eventdate<="2008-03-31" ~ "0708",
                          eventdate>="2008-04-01" & eventdate<="2009-03-31" ~ "0809",
                          eventdate>="2009-04-01" & eventdate<="2010-03-31" ~ "0910",
                          eventdate>="2010-04-01" & eventdate<="2011-03-31" ~ "1011",
                          eventdate>="2011-04-01" & eventdate<="2012-03-31" ~ "1112",
                          eventdate>="2012-04-01" & eventdate<="2013-03-31" ~ "1213",
                          eventdate>="2013-04-01" & eventdate<="2014-03-31" ~ "1314",
                          eventdate>="2014-04-01" & eventdate<="2015-03-31" ~ "1415",
                          eventdate>="2015-04-01" & eventdate<="2016-03-31" ~ "1516",
                          eventdate>="2016-04-01" & eventdate<="2017-03-31" ~ "1617",
                          eventdate>="2017-04-01" & eventdate<="2018-03-31" ~ "1718",
                          TRUE ~ "HELP"))
table(Ex$Year)

#Merge dates

Ex<-merge(x=Ex, y=CPRDActive, by="patid", all.x=TRUE, all.y=FALSE)

#If not active that year then drop
Ex<-Ex%>%
  mutate(drop=case_when(Year=="0405" & Act0405==0 ~ "drop",
                        Year=="0506" & Act0506==0 ~ "drop",
                        Year=="0607" & Act0607==0 ~ "drop",
                        Year=="0708" & Act0708==0 ~ "drop",
                        Year=="0809" & Act0809==0 ~ "drop",
                        Year=="0910" & Act0910==0 ~ "drop",
                        Year=="1011" & Act1011==0 ~ "drop",
                        Year=="1112" & Act1112==0 ~ "drop",
                        Year=="1213" & Act1213==0 ~ "drop",
                        Year=="1314" & Act1314==0 ~ "drop",
                        Year=="1415" & Act1415==0 ~ "drop",
                        Year=="1516" & Act1516==0 ~ "drop",
                        Year=="1617" & Act1617==0 ~ "drop",
                        Year=="1718" & Act1718==0 ~ "drop",
                        TRUE ~ "Keep"))%>%
  subset(drop=="Keep")

####Ever QOF####
QOF<-subset(Ex, QOF==1)
length(unique(QOF$patid))
NotQOF<-subset(Ex, QOF==0)
length(unique(NotQOF$patid))
           
####Table of different exceptions, ever reported####
table(Ex$Type)
table(Ex$Group)

ExSum<-Ex%>%
  group_by(patid)%>%
  select(patid, Group, Type)%>%
  distinct()%>%
  mutate(n=1)%>%
  pivot_wider(names_from=c("Group", "Type"), values_from="n", values_fill=0)

Group<-Ex%>%
  group_by(patid)%>%
  select(patid, Group)%>%
  distinct()%>%
  mutate(n=1)%>%
  pivot_wider(names_from=c("Group"), values_from="n", values_fill=0)

Type<-Ex%>%
  group_by(patid)%>%
  select(patid, Type)%>%
  distinct()%>%
  mutate(n=1)%>%
  pivot_wider(names_from=c("Type"), values_from="n", values_fill=0)

TypeMH<-Ex%>%
  group_by(patid)%>%
  subset(Group!="Other")%>%
  select(patid, Type)%>%
  distinct()%>%
  mutate(n=1)%>%
  pivot_wider(names_from=c("Type"), values_from="n", values_fill=0)%>%
  select(patid, MHdec=Declined, MHUns=Unsuitable, MHGen=` General`)  
  
###Got to here///// 

Patients<-select(CPRD, patid)
Patients<-merge(x=Patients, y=ExSum, by="patid", all.x=TRUE, all.y=TRUE)
Patients<-merge(x=Patients, y=Group, by="patid", all.x=TRUE, all.y=TRUE)
Patients<-merge(x=Patients, y=Type, by="patid", all.x=TRUE, all.y=TRUE)
Patients<-merge(x=Patients, y=TypeMH, by="patid", all.x=TRUE, all.y=TRUE)

Patients[is.na(Patients)] <- 0

#Who has both MH and other

#Table 1 - Cohort basics
MyVars<-names(Patients[,c(2:31)])

Table1<-CreateTableOne(vars=MyVars,  factorVars=MyVars, data=Patients, includeNA = TRUE)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=3)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=3, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Excepted.csv")

length(which(Patients$MH==1 & Patients$Other==1))

length(which(Patients$MH==1))
length(which(Patients$MH==1 & (Patients$Smoke==1|Patients$BMI==1|Patients$Alchol==1|Patients$BP==1|Patients$Cholesterol==1|Patients$Glucose==1)))
length(which(Patients$MH==1 | (Patients$Smoke==1|Patients$BMI==1|Patients$Alchol==1|Patients$BP==1|Patients$Cholesterol==1|Patients$Glucose==1)))

####Median number of MH exception reports####
Median<-Ex%>%
  subset(Group=="MH")%>%
  group_by(patid)%>%
  summarise(count=n())

summary(Median$count)

####Exception reporting in each year####

YearAll<-Ex%>%
  subset(Group=="MH")%>%
  select(patid, Year)%>%
  group_by(patid, Year)%>%
  mutate(n=n())%>%
  ungroup()%>%
  distinct()%>%
  group_by(patid)%>%
  arrange(Year)%>%
  pivot_wider(names_from = "Year", values_from = "n", values_fill= 0)

Exempt<-merge(x=YearAll, y=CPRDActive, by="patid", all.x=TRUE, all.y=TRUE)

Exempt[is.na(Exempt)] <- 0
Exempt[,2:15][Exempt[,2:15]>1]<-1

Years<-Exempt%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Active04=sum(c_across(16:29)),#Count up the number active years
         Ex04=sum(c_across(2:15)))%>%
  mutate(ExComplete04=case_when(Ex04==0 ~ "Never",
                                Active04==Ex04 ~ "Complete",
                                Active04>Ex04 ~ "Mixed",
                                TRUE ~ "Check"))
table(Years$ExComplete04)

Years$Prop<-Years$Ex04/Years$Active04*100
ExceptedOnly<-subset(Years, Ex04>0)
summary(ExceptedOnly$Prop)
length(which(ExceptedOnly$Prop>=50))

####Patient characteristics whole cohort####
CPRD<-merge(x=CPRD, y=Patients, by="patid", all.x=TRUE, all.y=TRUE)

MyVars<-c("source","AgeAtSMI", "AgeAtStart","gender", "ethn_missing", "LastDiag","country", "EnglandRegion", "FullIMD", "died", "AgeAtDeath",  "OtherQOF", "FU", "AnyPres", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_Exceptions_wholecohort.csv")

####Patient characteristics of those ever exception reported for mental health vs. never####
CPRD<-merge(x=CPRD, y=Patients, by="patid", all.x=TRUE, all.y=TRUE)

MyVars<-c("source","AgeAtSMI", "AgeAtStart","gender", "ethn_missing", "LastDiag","country", "EnglandRegion", "FullIMD", "died", "AgeAtDeath",  "OtherQOF", "FU", "AnyPres", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="MH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_Exceptions.csv")

####For those resident in England####
Eng<-subset(CPRD, country=="England")

MyVars<-c("EnglandRegion")
Table1<-CreateTableOne(vars=MyVars,  data=Eng, strata="MH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsEngland.csv")

IMD<-subset(CPRD, !is.na(FullIMD))

MyVars<-c("FullIMD")
Table1<-CreateTableOne(vars=MyVars,  data=IMD, strata="MH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsIMD.csv")

####Patient characteristics of those ever declined vs. never declined####
CPRD$EverDeclinedMH<-CPRD$MH_Declined

MyVars<-c("source","AgeAtSMI", "AgeAtStart","gender", "ethn_missing", "LastDiag","country", "EnglandRegion", "FullIMD", "died", "AgeAtDeath",  "OtherQOF", "FU", "AnyPres", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="EverDeclinedMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsDeclined.csv")

####For those resident in England####
Eng<-subset(CPRD, country=="England")

MyVars<-c("EnglandRegion")
Table1<-CreateTableOne(vars=MyVars,  data=Eng, strata="EverDeclinedMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsEnglandDeclined.csv")

IMD<-subset(CPRD, !is.na(FullIMD))

MyVars<-c("FullIMD")
Table1<-CreateTableOne(vars=MyVars,  data=IMD, strata="EverDeclinedMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsIMDDeclined.csv")

####Patient characteristics of those ever Unsuitable vs. never Unsuitable####
CPRD$EverUnsuitableMH<-CPRD$MH_Unsuitable

MyVars<-c("source","AgeAtSMI", "AgeAtStart","gender", "ethn_missing", "LastDiag","country", "EnglandRegion", "FullIMD", "died", "AgeAtDeath",  "OtherQOF", "FU", "AnyPres", "AnyScreen", "AllScreen")
Table1<-CreateTableOne(vars=MyVars,  data=CPRD, strata="EverUnsuitableMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)


write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsUnsuitable.csv")

####For those resident in England####
Eng<-subset(CPRD, country=="England")

MyVars<-c("EnglandRegion")
Table1<-CreateTableOne(vars=MyVars,  data=Eng, strata="EverUnsuitableMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsEnglandUnsuitable.csv")

IMD<-subset(CPRD, !is.na(FullIMD))

MyVars<-c("FullIMD")
Table1<-CreateTableOne(vars=MyVars,  data=IMD, strata="EverUnsuitableMH", includeNA = TRUE)
summary(Table1)
print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2)

Table1Exp <- print(Table1, printToggle = FALSE, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, catDigits=2, includeNA = TRUE)

write.csv(Table1Exp, file = "DatamindCPRD/Outputs/Table1_ExceptionsIMDUnsuitable.csv")

length(which(CPRD$EverUnsuitableMH==1 & CPRD$EverDeclinedMH==1))

####Split into different types of MH####

#Unique MH per year
Unique<-select(Ex, c(, 1, 7:9))
Unique<-subset(Unique, Group=="MH")
Unique<-distinct(Unique)

Unique<-merge(x=CPRDActive, y=Unique, by="patid", all.x=TRUE, all.y=FALSE)

Unsuitable<-subset(Unique, Type =="Unsuitable")
Declined<-subset(Unique, Type =="Declined")
Any<-select(Unique, -Type)
Any<-distinct(Any)

Unique<-select(Unique, -Type, -Group, -Year)
Unique<-distinct(Unique)

AnyEx<-Any%>%
  group_by(patid, Year)%>%
  subset(!is.na(Year))%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)%>%
  mutate(Type="Any")

Unsuitable<-Unsuitable%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)%>%
  mutate(Type="Unsuitable")

Declined<-Declined%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)%>%
  mutate(Type="Declined")

All<-rbind(AnyEx, Unsuitable, Declined)

####Main graph####
AllAny<-All%>%
  group_by(Type)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Unique%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")

Active<-rbind(AllActive, AllActive, AllActive)
AllType<-cbind(Active, AllAny)

AllType$Year<-str_replace(AllType$ActYear, "Act", "")

AllType<-select(AllType, strat=Type, Year, AnyCount, ActCount)

TypeCI<-BinomCI(AllType$AnyCount, AllType$ActCount)*100

AllType<-cbind(AllType, TypeCI)

AllType$Year<-case_when(AllType$Year=="0405" ~ "2004-2005",
                        AllType$Year=="0506" ~ "2005-2006",
                        AllType$Year=="0607" ~ "2006-2007",
                        AllType$Year=="0708" ~ "2007-2008",
                        AllType$Year=="0809" ~ "2008-2009",
                        AllType$Year=="0910" ~ "2009-2010",
                        AllType$Year=="1011" ~ "2010-2011",
                        AllType$Year=="1112" ~ "2011-2012",
                        AllType$Year=="1213" ~ "2012-2013",
                        AllType$Year=="1314" ~ "2013-2014",
                        AllType$Year=="1415" ~ "2014-2015",
                        AllType$Year=="1516" ~ "2015-2016",
                        AllType$Year=="1617" ~ "2016-2017",
                        AllType$Year=="1718" ~ "2017-2018")

AllType$strat[AllType$strat=="Declined"]<-"Informed dissent"


ggplot()+
  geom_line(data = AllType, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllType, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.15, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 16), breaks = seq(0,16,by = 2))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="Exception type")

ggsave("DatamindCPRD/Outputs/ExceptionAll.PNG", width = 10,   height = 10)
ggsave("DatamindCPRD/Outputs/ExceptionAll.pdf", width = 10,   height = 10)

####Stratification - Any exception report####

#Create new vars
NewVars<-select(CPRD, patid, country, LastDiag, gender, yob, OtherQOF, pracid, ethn_missing, FullIMD)
Unique<-merge(x=Unique, y=NewVars, by="patid", all.x=TRUE, all.y=TRUE)

AnyEx<-merge(x=AnyEx, y=Unique, by="patid", all.x=TRUE, all.y=TRUE)
AnyEx[,2:15][is.na(AnyEx[,2:15])]<-0


#By country
#Calculate totals

AllAny<-AnyEx%>%
  group_by(country)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AnyEx%>%
  group_by(country)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-country)

AllCountry<-cbind(AllActive, AllAny)

AllCountry$Year<-str_replace(AllCountry$ActYear, "Act", "")

AllCountry<-select(AllCountry, strat=country, Year, AnyCount, ActCount)
AllCountry$Group<-"Country"

CountryCI<-BinomCI(AllCountry$AnyCount, AllCountry$ActCount)*100

AllCountry<-cbind(AllCountry, CountryCI)

AllCountry$Year<-case_when(AllCountry$Year=="0405" ~ "2004-2005",
                           AllCountry$Year=="0506" ~ "2005-2006",
                           AllCountry$Year=="0607" ~ "2006-2007",
                           AllCountry$Year=="0708" ~ "2007-2008",
                           AllCountry$Year=="0809" ~ "2008-2009",
                           AllCountry$Year=="0910" ~ "2009-2010",
                           AllCountry$Year=="1011" ~ "2010-2011",
                           AllCountry$Year=="1112" ~ "2011-2012",
                           AllCountry$Year=="1213" ~ "2012-2013",
                           AllCountry$Year=="1314" ~ "2013-2014",
                           AllCountry$Year=="1415" ~ "2014-2015",
                           AllCountry$Year=="1516" ~ "2015-2016",
                           AllCountry$Year=="1617" ~ "2016-2017",
                           AllCountry$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllCountry, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllCountry, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="Country")

ggsave("DatamindCPRD/Outputs/ExceptionCountry.PNG", width = 10,   height = 10)

#By gender
AllAny<-AnyEx%>%
  group_by(gender)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AnyEx%>%
  group_by(gender)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-gender)

AllGender<-cbind(AllActive, AllAny)

AllGender$Year<-str_replace(AllGender$ActYear, "Act", "")

AllGender<-select(AllGender, strat=gender, Year, AnyCount, ActCount)
AllGender$Group<-"Sex"

GenderCI<-BinomCI(AllGender$AnyCount, AllGender$ActCount)*100

AllGender<-cbind(AllGender, GenderCI)

AllGender$Year<-case_when(AllGender$Year=="0405" ~ "2004-2005",
                          AllGender$Year=="0506" ~ "2005-2006",
                          AllGender$Year=="0607" ~ "2006-2007",
                          AllGender$Year=="0708" ~ "2007-2008",
                          AllGender$Year=="0809" ~ "2008-2009",
                          AllGender$Year=="0910" ~ "2009-2010",
                          AllGender$Year=="1011" ~ "2010-2011",
                          AllGender$Year=="1112" ~ "2011-2012",
                          AllGender$Year=="1213" ~ "2012-2013",
                          AllGender$Year=="1314" ~ "2013-2014",
                          AllGender$Year=="1415" ~ "2014-2015",
                          AllGender$Year=="1516" ~ "2015-2016",
                          AllGender$Year=="1617" ~ "2016-2017",
                          AllGender$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllGender, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllGender, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="Sex")

ggsave("DatamindCPRD/Outputs/ExceptionSex.PNG", width = 10,   height = 10)

#By ethnicity
AllAny<-AnyEx%>%
  group_by(ethn_missing)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AnyEx%>%
  group_by(ethn_missing)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-ethn_missing)

AllEthn<-cbind(AllActive, AllAny)

AllEthn$Year<-str_replace(AllEthn$ActYear, "Act", "")

AllEthn<-select(AllEthn, strat=ethn_missing, Year, AnyCount, ActCount)
AllEthn$Group<-"Ethnicity"

EthnCI<-BinomCI(AllEthn$AnyCount, AllEthn$ActCount)*100

AllEthn<-cbind(AllEthn, EthnCI)

AllEthn$Year<-case_when(AllEthn$Year=="0405" ~ "2004-2005",
                        AllEthn$Year=="0506" ~ "2005-2006",
                        AllEthn$Year=="0607" ~ "2006-2007",
                        AllEthn$Year=="0708" ~ "2007-2008",
                        AllEthn$Year=="0809" ~ "2008-2009",
                        AllEthn$Year=="0910" ~ "2009-2010",
                        AllEthn$Year=="1011" ~ "2010-2011",
                        AllEthn$Year=="1112" ~ "2011-2012",
                        AllEthn$Year=="1213" ~ "2012-2013",
                        AllEthn$Year=="1314" ~ "2013-2014",
                        AllEthn$Year=="1415" ~ "2014-2015",
                        AllEthn$Year=="1516" ~ "2015-2016",
                        AllEthn$Year=="1617" ~ "2016-2017",
                        AllEthn$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllEthn, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="Ethnicity")

ggsave("DatamindCPRD/Outputs/ExceptionEthnicity.PNG", width = 10,   height = 10)

#Without CI
ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Ethnicity")

####Now think about time varying ones####
####Age####

Age<-AnyEx%>%
  mutate(Age0405 = case_when(2004-yob>=18 & 2004-yob<=30~ "18-30",
                             2004-yob>30 & 2004-yob<=40~ "31-40",
                             2004-yob>40 & 2004-yob<=50~ "41-50",
                             2004-yob>50 & 2004-yob<=60~ "51-60",
                             2004-yob>60 & 2004-yob<=70~ "61-70",
                             2004-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0506 = case_when(2005-yob>=18 & 2005-yob<=30~ "18-30",
                             2005-yob>30 & 2005-yob<=40~ "31-40",
                             2005-yob>40 & 2005-yob<=50~ "41-50",
                             2005-yob>50 & 2005-yob<=60~ "51-60",
                             2005-yob>60 & 2005-yob<=70~ "61-70",
                             2005-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0607 = case_when(2006-yob>=18 & 2006-yob<=30~ "18-30",
                             2006-yob>30 & 2006-yob<=40~ "31-40",
                             2006-yob>40 & 2006-yob<=50~ "41-50",
                             2006-yob>50 & 2006-yob<=60~ "51-60",
                             2006-yob>60 & 2006-yob<=70~ "61-70",
                             2006-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0708 = case_when(2007-yob>=18 & 2007-yob<=30~ "18-30",
                             2007-yob>30 & 2007-yob<=40~ "31-40",
                             2007-yob>40 & 2007-yob<=50~ "41-50",
                             2007-yob>50 & 2007-yob<=60~ "51-60",
                             2007-yob>60 & 2007-yob<=70~ "61-70",
                             2007-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0809 = case_when(2008-yob>=18 & 2008-yob<=30~ "18-30",
                             2008-yob>30 & 2008-yob<=40~ "31-40",
                             2008-yob>40 & 2008-yob<=50~ "41-50",
                             2008-yob>50 & 2008-yob<=60~ "51-60",
                             2008-yob>60 & 2008-yob<=70~ "61-70",
                             2008-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0910 = case_when(2009-yob>=18 & 2009-yob<=30~ "18-30",
                             2009-yob>30 & 2009-yob<=40~ "31-40",
                             2009-yob>40 & 2009-yob<=50~ "41-50",
                             2009-yob>50 & 2009-yob<=60~ "51-60",
                             2009-yob>60 & 2009-yob<=70~ "61-70",
                             2009-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1011 = case_when(2010-yob>=18 & 2010-yob<=30~ "18-30",
                             2010-yob>30 & 2010-yob<=40~ "31-40",
                             2010-yob>40 & 2010-yob<=50~ "41-50",
                             2010-yob>50 & 2010-yob<=60~ "51-60",
                             2010-yob>60 & 2010-yob<=70~ "61-70",
                             2010-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1112 = case_when(2011-yob>=18 & 2011-yob<=30~ "18-30",
                             2011-yob>30 & 2011-yob<=40~ "31-40",
                             2011-yob>40 & 2011-yob<=50~ "41-50",
                             2011-yob>50 & 2011-yob<=60~ "51-60",
                             2011-yob>60 & 2011-yob<=70~ "61-70",
                             2011-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1213 = case_when(2012-yob>=18 & 2012-yob<=30~ "18-30",
                             2012-yob>30 & 2012-yob<=40~ "31-40",
                             2012-yob>40 & 2012-yob<=50~ "41-50",
                             2012-yob>50 & 2012-yob<=60~ "51-60",
                             2012-yob>60 & 2012-yob<=70~ "61-70",
                             2012-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1314 = case_when(2013-yob>=18 & 2013-yob<=30~ "18-30",
                             2013-yob>30 & 2013-yob<=40~ "31-40",
                             2013-yob>40 & 2013-yob<=50~ "41-50",
                             2013-yob>50 & 2013-yob<=60~ "51-60",
                             2013-yob>60 & 2013-yob<=70~ "61-70",
                             2013-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1415 = case_when(2014-yob>=18 & 2014-yob<=30~ "18-30",
                             2014-yob>30 & 2014-yob<=40~ "31-40",
                             2014-yob>40 & 2014-yob<=50~ "41-50",
                             2014-yob>50 & 2014-yob<=60~ "51-60",
                             2014-yob>60 & 2014-yob<=70~ "61-70",
                             2014-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1516 = case_when(2015-yob>=18 & 2015-yob<=30~ "18-30",
                             2015-yob>30 & 2015-yob<=40~ "31-40",
                             2015-yob>40 & 2015-yob<=50~ "41-50",
                             2015-yob>50 & 2015-yob<=60~ "51-60",
                             2015-yob>60 & 2015-yob<=70~ "61-70",
                             2015-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1617 = case_when(2016-yob>=18 & 2016-yob<=30~ "18-30",
                             2016-yob>30 & 2016-yob<=40~ "31-40",
                             2016-yob>40 & 2016-yob<=50~ "41-50",
                             2016-yob>50 & 2016-yob<=60~ "51-60",
                             2016-yob>60 & 2016-yob<=70~ "61-70",
                             2016-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1718 = case_when(2017-yob>=18 & 2017-yob<=30~ "18-30",
                             2017-yob>30 & 2017-yob<=40~ "31-40",
                             2017-yob>40 & 2017-yob<=50~ "41-50",
                             2017-yob>50 & 2017-yob<=60~ "51-60",
                             2017-yob>60 & 2017-yob<=70~ "61-70",
                             2017-yob>70 ~ ">70",
                             TRUE ~ "HELP"))


#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAge<-subset(Age, Act0405==1)%>%
  select(patid, Age=Age0405, Any=Any0405)%>%
  group_by(Age)%>%
  summarise(ActCount=n(),AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
NewFile<-subset(Age, Age[paste0("Act", YearVar[i])]==1)
NewFile<-select(NewFile, patid, Age=paste0("Age", YearVar[i]), Any=paste0("Any", YearVar[i]))

NewFile<-NewFile%>%
    group_by(Age)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
AllAge<-rbind(AllAge, NewFile)
} 

AllAge$AgeF<-factor(AllAge$Age, levels=c("18-30", "31-40", "41-50", "51-60", "61-70", ">70"))

AllAge<-select(AllAge, strat=AgeF, Year, AnyCount, ActCount)
AllAge$Group<-"Age"

AllAgeCI<-BinomCI(AllAge$AnyCount, AllAge$ActCount)*100

AllAge<-cbind(AllAge, AllAgeCI)

AllAge$Year<-case_when(AllAge$Year=="0405" ~ "2004-2005",
                       AllAge$Year=="0506" ~ "2005-2006",
                       AllAge$Year=="0607" ~ "2006-2007",
                       AllAge$Year=="0708" ~ "2007-2008",
                       AllAge$Year=="0809" ~ "2008-2009",
                       AllAge$Year=="0910" ~ "2009-2010",
                       AllAge$Year=="1011" ~ "2010-2011",
                       AllAge$Year=="1112" ~ "2011-2012",
                       AllAge$Year=="1213" ~ "2012-2013",
                       AllAge$Year=="1314" ~ "2013-2014",
                       AllAge$Year=="1415" ~ "2014-2015",
                       AllAge$Year=="1516" ~ "2015-2016",
                       AllAge$Year=="1617" ~ "2016-2017",
                       AllAge$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllAge, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAge, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="Age")

ggsave("DatamindCPRD/Outputs/ExceptionAge.PNG", width = 10,   height = 10)

####Diagnosis long####
load("VariableExtracts/CPRD2018/MergedObs/SMIObs.Rdata")
SMIObs<-subset(SMIObs, patid %in% CPRD$patid)

Diag<-SMIObs%>%
  group_by(patid)%>%
  mutate(DiagYear = case_when(SMIDate<="2005-03-31" ~ "0405",
                            SMIDate>="2005-04-01" & SMIDate<="2006-03-31" ~ "0506",
                            SMIDate>="2006-04-01" & SMIDate<="2007-03-31" ~ "0607",
                            SMIDate>="2007-04-01" & SMIDate<="2008-03-31" ~ "0708",
                            SMIDate>="2008-04-01" & SMIDate<="2009-03-31" ~ "0809",
                            SMIDate>="2009-04-01" & SMIDate<="2010-03-31" ~ "0910",
                            SMIDate>="2010-04-01" & SMIDate<="2011-03-31" ~ "1011",
                            SMIDate>="2011-04-01" & SMIDate<="2012-03-31" ~ "1112",
                            SMIDate>="2012-04-01" & SMIDate<="2013-03-31" ~ "1213",
                            SMIDate>="2013-04-01" & SMIDate<="2014-03-31" ~ "1314",
                            SMIDate>="2014-04-01" & SMIDate<="2015-03-31" ~ "1415",
                            SMIDate>="2015-04-01" & SMIDate<="2016-03-31" ~ "1516",
                            SMIDate>="2016-04-01" & SMIDate<="2017-03-31" ~ "1617",
                            SMIDate>="2017-04-01" & SMIDate<="2018-03-31" ~ "1718",
                            TRUE ~ NA))
Diag1<-Diag%>%
  mutate(Order=case_when(diag=="schizophrenia" ~ 1,
                         diag=="bipolar" ~ 2,
                         diag=="other psychosis" ~ 3,
                         TRUE ~ 0))%>%
  group_by(patid, DiagYear)%>%
  mutate(Priority=min(Order))%>%
  filter(Priority==Order & !is.na(DiagYear))%>%
  select(patid, DiagYear, diag)%>%
  distinct()%>%
  arrange(DiagYear)%>%
  pivot_wider(names_from=DiagYear, names_prefix = "Diag",  values_from = diag)

Diag2<-Diag1%>%
  pivot_longer(cols= c(2:15), names_to = "DiagYear", values_to = "diag")%>%
  group_by(patid)%>%
  mutate(diag=case_when(is.na(diag) ~ lag(diag),
                        TRUE ~ diag))%>%
  fill(diag, .direction=c("down"))

Diag3<-Diag2%>%
  subset(DiagYear!="Diag0000")%>%
  pivot_wider(names_from=DiagYear, values_from = diag)

Diag<-merge(x=AnyEx, y=Diag3, by="patid", all.x=TRUE, All.y=FALSE)

DiagVar<-vars_select(names(Diag), starts_with('Diag', ignore.case = TRUE))

#Run it by year
AllDiag<-subset(Diag, Act0405==1)%>%
  select(patid, Diag=Diag0405, Any=Any0405)%>%
  group_by(Diag)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(DiagVar))) {
  NewFile<-subset(Diag, Diag[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Diag=paste0("Diag", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Diag)%>%
    summarise(ActCount=n(), AnyCount=sum(Any),Year=as.character(paste0(YearVar[i])))  
  AllDiag<-rbind(AllDiag, NewFile)
} 

AllDiag$Diag<-as.factor(AllDiag$Diag)

levels(AllDiag$Diag)<-c("Bipolar disorder", "Other psychoses", "Schizophrenia")
AllDiag$Diag <- factor(AllDiag$Diag, levels = c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

AllDiag<-select(AllDiag, strat=Diag, Year, AnyCount, ActCount)
AllDiag$Group<-"DiagTV"

AllDiagCI<-BinomCI(AllDiag$AnyCount, AllDiag$ActCount)*100

AllDiag<-cbind(AllDiag, AllDiagCI)

AllDiag$Year<-case_when(AllDiag$Year=="0405" ~ "2004-2005",
                        AllDiag$Year=="0506" ~ "2005-2006",
                        AllDiag$Year=="0607" ~ "2006-2007",
                        AllDiag$Year=="0708" ~ "2007-2008",
                        AllDiag$Year=="0809" ~ "2008-2009",
                        AllDiag$Year=="0910" ~ "2009-2010",
                        AllDiag$Year=="1011" ~ "2010-2011",
                        AllDiag$Year=="1112" ~ "2011-2012",
                        AllDiag$Year=="1213" ~ "2012-2013",
                        AllDiag$Year=="1314" ~ "2013-2014",
                       AllDiag$Year=="1415" ~ "2014-2015",
                       AllDiag$Year=="1516" ~ "2015-2016",
                       AllDiag$Year=="1617" ~ "2016-2017",
                       AllDiag$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllDiag, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllDiag, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="SMI diagnosis")

ggsave("DatamindCPRD/Outputs/ExceptionDiag.PNG", width = 10,   height = 10)

####Antipsychotics####
CPRD$FirstAnyPres<-pmin(CPRD$FirstAPDate, CPRD$FirstBPDate, CPRD$FirstEvidenceAPDate, CPRD$FirstEvidenceBPDate, na.rm=TRUE)
APFields<-select(CPRD, patid, FirstAnyPres)
AnyEx<-merge(x=AnyEx, y=APFields, by="patid", all.x=TRUE, all.y=FALSE)

AP<-AnyEx%>%
  mutate(AP0405 = case_when(FirstAnyPres<="2005-03-31" ~ 1,
                               TRUE ~ 0),
         AP0506 = case_when(FirstAnyPres<="2006-03-31" ~ 1,
                               TRUE ~ 0),
         AP0607 = case_when(FirstAnyPres<="2007-03-31" ~ 1,
                               TRUE ~ 0),
         AP0708 = case_when(FirstAnyPres<="2008-03-31" ~ 1,
                               TRUE ~ 0), 
         AP0809 = case_when(FirstAnyPres<="2009-03-31"  ~ 1,
                               TRUE ~ 0),
         AP0910 = case_when(FirstAnyPres<="2010-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1011 = case_when(FirstAnyPres<="2011-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1112 = case_when(FirstAnyPres<="2012-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1213 = case_when(FirstAnyPres<="2013-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1314 = case_when(FirstAnyPres<="2014-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1415 = case_when(FirstAnyPres<="2015-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1516 = case_when(FirstAnyPres<="2016-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1617 = case_when(FirstAnyPres<="2017-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1718 = case_when(FirstAnyPres<="2018-03-31" ~ 1,
                               TRUE ~ 0))

#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAP<-subset(AP, Act0405==1)%>%
  select(patid, AP=AP0405, Any=Any0405)%>%
  group_by(AP)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
  NewFile<-subset(AP, AP[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, AP=paste0("AP", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(AP)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllAP<-rbind(AllAP, NewFile)
} 

AllAP$AP<-as.factor(AllAP$AP)
levels(AllAP$AP)<-c("No", "Yes")

AllAP<-select(AllAP, strat=AP, Year, AnyCount, ActCount)
AllAP$Group<-"APTV"

AllAPCI<-BinomCI(AllAP$AnyCount, AllAP$ActCount)*100

AllAP<-cbind(AllAP, AllAPCI)

AllAP$Year<-case_when(AllAP$Year=="0405" ~ "2004-2005",
                      AllAP$Year=="0506" ~ "2005-2006",
                      AllAP$Year=="0607" ~ "2006-2007",
                      AllAP$Year=="0708" ~ "2007-2008",
                      AllAP$Year=="0809" ~ "2008-2009",
                      AllAP$Year=="0910" ~ "2009-2010",
                      AllAP$Year=="1011" ~ "2010-2011",
                      AllAP$Year=="1112" ~ "2011-2012",
                      AllAP$Year=="1213" ~ "2012-2013",
                      AllAP$Year=="1314" ~ "2013-2014",
                      AllAP$Year=="1415" ~ "2014-2015",
                      AllAP$Year=="1516" ~ "2015-2016",
                      AllAP$Year=="1617" ~ "2016-2017",
                      AllAP$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllAP, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAP, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="on antipsychotics/mood stabilisers")

ggsave("DatamindCPRD/Outputs/ExceptionAP.PNG", width = 10,   height = 10)

####For those with deprivation data####
AllAny<-AnyEx%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-AnyEx%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-FullIMD)

AllIMD<-cbind(AllActive, AllAny)

AllIMD$Year<-str_replace(AllIMD$ActYear, "Act", "")

AllIMD<-select(AllIMD, strat=FullIMD, Year, AnyCount, ActCount)
AllIMD$Group<-"IMD"

IMDCI<-BinomCI(AllIMD$AnyCount, AllIMD$ActCount)*100

AllIMD<-cbind(AllIMD, IMDCI)

AllIMD$Year<-case_when(AllIMD$Year=="0405" ~ "2004-2005",
                       AllIMD$Year=="0506" ~ "2005-2006",
                       AllIMD$Year=="0607" ~ "2006-2007",
                       AllIMD$Year=="0708" ~ "2007-2008",
                       AllIMD$Year=="0809" ~ "2008-2009",
                       AllIMD$Year=="0910" ~ "2009-2010",
                       AllIMD$Year=="1011" ~ "2010-2011",
                       AllIMD$Year=="1112" ~ "2011-2012",
                       AllIMD$Year=="1213" ~ "2012-2013",
                       AllIMD$Year=="1314" ~ "2013-2014",
                       AllIMD$Year=="1415" ~ "2014-2015",
                       AllIMD$Year=="1516" ~ "2015-2016",
                       AllIMD$Year=="1617" ~ "2016-2017",
                       AllIMD$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllIMD, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllIMD, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", fill="IMD")

ggsave("DatamindCPRD/Outputs/ExceptionIMD.PNG", width = 10,   height = 10)

All<-rbind(AllCountry, AllGender, AllEthn, AllAge, AllDiag, AllAP, AllIMD)
write.csv(All, file="DatamindCPRD/Outputs/ExceptAll.csv")

#####Stratified for declined####
#Create new vars
Declined<-merge(x=Declined, y=Unique, by="patid", all.x=TRUE, all.y=TRUE)
Declined[,2:15][is.na(Declined[,2:15])]<-0

#By country
#Calculate totals

AllAny<-Declined%>%
  group_by(country)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Declined%>%
  group_by(country)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-country)

AllCountry<-cbind(AllActive, AllAny)

AllCountry$Year<-str_replace(AllCountry$ActYear, "Act", "")

AllCountry<-select(AllCountry, strat=country, Year, AnyCount, ActCount)
AllCountry$Group<-"Country"

CountryCI<-BinomCI(AllCountry$AnyCount, AllCountry$ActCount)*100

AllCountry<-cbind(AllCountry, CountryCI)

ggplot()+
  geom_line(data = AllCountry, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllCountry, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Country")

#By gender
AllAny<-Declined%>%
  group_by(gender)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Declined%>%
  group_by(gender)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-gender)

AllGender<-cbind(AllActive, AllAny)

AllGender$Year<-str_replace(AllGender$ActYear, "Act", "")

AllGender<-select(AllGender, strat=gender, Year, AnyCount, ActCount)
AllGender$Group<-"Sex"

GenderCI<-BinomCI(AllGender$AnyCount, AllGender$ActCount)*100

AllGender<-cbind(AllGender, GenderCI)

ggplot()+
  geom_line(data = AllGender, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllGender, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Sex")

#By ethnicity
AllAny<-Declined%>%
  group_by(ethn_missing)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Declined%>%
  group_by(ethn_missing)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-ethn_missing)

AllEthn<-cbind(AllActive, AllAny)

AllEthn$Year<-str_replace(AllEthn$ActYear, "Act", "")

AllEthn<-select(AllEthn, strat=ethn_missing, Year, AnyCount, ActCount)
AllEthn$Group<-"Ethnicity"

EthnCI<-BinomCI(AllEthn$AnyCount, AllEthn$ActCount)*100

AllEthn<-cbind(AllEthn, EthnCI)

ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllEthn, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Ethnicity")

#Without CI
ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Ethnicity")

####Now think about time varying ones####
####Age####

Age<-Declined%>%
  mutate(Age0405 = case_when(2004-yob>=18 & 2004-yob<=30~ "18-30",
                             2004-yob>30 & 2004-yob<=40~ "31-40",
                             2004-yob>40 & 2004-yob<=50~ "41-50",
                             2004-yob>50 & 2004-yob<=60~ "51-60",
                             2004-yob>60 & 2004-yob<=70~ "61-70",
                             2004-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0506 = case_when(2005-yob>=18 & 2005-yob<=30~ "18-30",
                             2005-yob>30 & 2005-yob<=40~ "31-40",
                             2005-yob>40 & 2005-yob<=50~ "41-50",
                             2005-yob>50 & 2005-yob<=60~ "51-60",
                             2005-yob>60 & 2005-yob<=70~ "61-70",
                             2005-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0607 = case_when(2006-yob>=18 & 2006-yob<=30~ "18-30",
                             2006-yob>30 & 2006-yob<=40~ "31-40",
                             2006-yob>40 & 2006-yob<=50~ "41-50",
                             2006-yob>50 & 2006-yob<=60~ "51-60",
                             2006-yob>60 & 2006-yob<=70~ "61-70",
                             2006-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0708 = case_when(2007-yob>=18 & 2007-yob<=30~ "18-30",
                             2007-yob>30 & 2007-yob<=40~ "31-40",
                             2007-yob>40 & 2007-yob<=50~ "41-50",
                             2007-yob>50 & 2007-yob<=60~ "51-60",
                             2007-yob>60 & 2007-yob<=70~ "61-70",
                             2007-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0809 = case_when(2008-yob>=18 & 2008-yob<=30~ "18-30",
                             2008-yob>30 & 2008-yob<=40~ "31-40",
                             2008-yob>40 & 2008-yob<=50~ "41-50",
                             2008-yob>50 & 2008-yob<=60~ "51-60",
                             2008-yob>60 & 2008-yob<=70~ "61-70",
                             2008-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0910 = case_when(2009-yob>=18 & 2009-yob<=30~ "18-30",
                             2009-yob>30 & 2009-yob<=40~ "31-40",
                             2009-yob>40 & 2009-yob<=50~ "41-50",
                             2009-yob>50 & 2009-yob<=60~ "51-60",
                             2009-yob>60 & 2009-yob<=70~ "61-70",
                             2009-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1011 = case_when(2010-yob>=18 & 2010-yob<=30~ "18-30",
                             2010-yob>30 & 2010-yob<=40~ "31-40",
                             2010-yob>40 & 2010-yob<=50~ "41-50",
                             2010-yob>50 & 2010-yob<=60~ "51-60",
                             2010-yob>60 & 2010-yob<=70~ "61-70",
                             2010-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1112 = case_when(2011-yob>=18 & 2011-yob<=30~ "18-30",
                             2011-yob>30 & 2011-yob<=40~ "31-40",
                             2011-yob>40 & 2011-yob<=50~ "41-50",
                             2011-yob>50 & 2011-yob<=60~ "51-60",
                             2011-yob>60 & 2011-yob<=70~ "61-70",
                             2011-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1213 = case_when(2012-yob>=18 & 2012-yob<=30~ "18-30",
                             2012-yob>30 & 2012-yob<=40~ "31-40",
                             2012-yob>40 & 2012-yob<=50~ "41-50",
                             2012-yob>50 & 2012-yob<=60~ "51-60",
                             2012-yob>60 & 2012-yob<=70~ "61-70",
                             2012-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1314 = case_when(2013-yob>=18 & 2013-yob<=30~ "18-30",
                             2013-yob>30 & 2013-yob<=40~ "31-40",
                             2013-yob>40 & 2013-yob<=50~ "41-50",
                             2013-yob>50 & 2013-yob<=60~ "51-60",
                             2013-yob>60 & 2013-yob<=70~ "61-70",
                             2013-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1415 = case_when(2014-yob>=18 & 2014-yob<=30~ "18-30",
                             2014-yob>30 & 2014-yob<=40~ "31-40",
                             2014-yob>40 & 2014-yob<=50~ "41-50",
                             2014-yob>50 & 2014-yob<=60~ "51-60",
                             2014-yob>60 & 2014-yob<=70~ "61-70",
                             2014-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1516 = case_when(2015-yob>=18 & 2015-yob<=30~ "18-30",
                             2015-yob>30 & 2015-yob<=40~ "31-40",
                             2015-yob>40 & 2015-yob<=50~ "41-50",
                             2015-yob>50 & 2015-yob<=60~ "51-60",
                             2015-yob>60 & 2015-yob<=70~ "61-70",
                             2015-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1617 = case_when(2016-yob>=18 & 2016-yob<=30~ "18-30",
                             2016-yob>30 & 2016-yob<=40~ "31-40",
                             2016-yob>40 & 2016-yob<=50~ "41-50",
                             2016-yob>50 & 2016-yob<=60~ "51-60",
                             2016-yob>60 & 2016-yob<=70~ "61-70",
                             2016-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1718 = case_when(2017-yob>=18 & 2017-yob<=30~ "18-30",
                             2017-yob>30 & 2017-yob<=40~ "31-40",
                             2017-yob>40 & 2017-yob<=50~ "41-50",
                             2017-yob>50 & 2017-yob<=60~ "51-60",
                             2017-yob>60 & 2017-yob<=70~ "61-70",
                             2017-yob>70 ~ ">70",
                             TRUE ~ "HELP"))


#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAge<-subset(Age, Act0405==1)%>%
  select(patid, Age=Age0405, Any=Any0405)%>%
  group_by(Age)%>%
  summarise(ActCount=n(),AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
NewFile<-subset(Age, Age[paste0("Act", YearVar[i])]==1)
NewFile<-select(NewFile, patid, Age=paste0("Age", YearVar[i]), Any=paste0("Any", YearVar[i]))

NewFile<-NewFile%>%
    group_by(Age)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
AllAge<-rbind(AllAge, NewFile)
} 

AllAge$AgeF<-as.factor(AllAge$Age)
levels(AllAge$Age)<-c("18-30", "31-40", "41-50", "51-60", "61-70", ">70")

AllAge<-select(AllAge, strat=Age, Year, AnyCount, ActCount)
AllAge$Group<-"Age"

AllAgeCI<-BinomCI(AllAge$AnyCount, AllAge$ActCount)*100

AllAge<-cbind(AllAge, AllAgeCI)

ggplot()+
  geom_line(data = AllAge, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAge, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Age")

####Diagnosis long####

Diag<-merge(x=Declined, y=Diag3, by="patid", all.x=TRUE, All.y=FALSE)

DiagVar<-vars_select(names(Diag), starts_with('Diag', ignore.case = TRUE))

#Run it by year
AllDiag<-subset(Diag, Act0405==1)%>%
  select(patid, Diag=Diag0405, Any=Any0405)%>%
  group_by(Diag)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(DiagVar))) {
  NewFile<-subset(Diag, Diag[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Diag=paste0("Diag", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Diag)%>%
    summarise(ActCount=n(), AnyCount=sum(Any),Year=as.character(paste0(YearVar[i])))  
  AllDiag<-rbind(AllDiag, NewFile)
} 

AllDiag$Diag<-as.factor(AllDiag$Diag)

levels(AllDiag$Diag)<-c("Bipolar disorder", "Other psychoses", "Schizophrenia")
AllDiag$Diag <- factor(AllDiag$Diag, levels = c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

AllDiag<-select(AllDiag, strat=Diag, Year, AnyCount, ActCount)
AllDiag$Group<-"DiagTV"

AllDiagCI<-BinomCI(AllDiag$AnyCount, AllDiag$ActCount)*100

AllDiag<-cbind(AllDiag, AllDiagCI)

ggplot()+
  geom_line(data = AllDiag, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllDiag, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="SMI diagnosis")

####Antipsychotics####
Declined<-merge(x=Declined, y=APFields, by="patid", all.x=TRUE, all.y=FALSE)

AP<-Declined%>%
  mutate(AP0405 = case_when(FirstAnyPres<="2005-03-31" ~ 1,
                               TRUE ~ 0),
         AP0506 = case_when(FirstAnyPres<="2006-03-31" ~ 1,
                               TRUE ~ 0),
         AP0607 = case_when(FirstAnyPres<="2007-03-31" ~ 1,
                               TRUE ~ 0),
         AP0708 = case_when(FirstAnyPres<="2008-03-31" ~ 1,
                               TRUE ~ 0), 
         AP0809 = case_when(FirstAnyPres<="2009-03-31"  ~ 1,
                               TRUE ~ 0),
         AP0910 = case_when(FirstAnyPres<="2010-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1011 = case_when(FirstAnyPres<="2011-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1112 = case_when(FirstAnyPres<="2012-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1213 = case_when(FirstAnyPres<="2013-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1314 = case_when(FirstAnyPres<="2014-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1415 = case_when(FirstAnyPres<="2015-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1516 = case_when(FirstAnyPres<="2016-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1617 = case_when(FirstAnyPres<="2017-03-31" ~ 1,
                               TRUE ~ 0), 
         AP1718 = case_when(FirstAnyPres<="2018-03-31" ~ 1,
                               TRUE ~ 0))

#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAP<-subset(AP, Act0405==1)%>%
  select(patid, AP=AP0405, Any=Any0405)%>%
  group_by(AP)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
  NewFile<-subset(AP, AP[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, AP=paste0("AP", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(AP)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllAP<-rbind(AllAP, NewFile)
} 

AllAP$AP<-as.factor(AllAP$AP)
levels(AllAP$AP)<-c("No", "Yes")

AllAP<-select(AllAP, strat=AP, Year, AnyCount, ActCount)
AllAP$Group<-"APTV"

AllAPCI<-BinomCI(AllAP$AnyCount, AllAP$ActCount)*100

AllAP<-cbind(AllAP, AllAPCI)

ggplot()+
  geom_line(data = AllAP, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAP, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.3, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0,26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="On antipsychotics/mood stabilisers")

####For those with deprivation data####
AllAny<-Declined%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Declined%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-FullIMD)

AllIMD<-cbind(AllActive, AllAny)

AllIMD$Year<-str_replace(AllIMD$ActYear, "Act", "")

AllIMD<-select(AllIMD, strat=FullIMD, Year, AnyCount, ActCount)
AllIMD$Group<-"IMD"

IMDCI<-BinomCI(AllIMD$AnyCount, AllIMD$ActCount)*100

AllIMD<-cbind(AllIMD, IMDCI)

ggplot()+
  geom_line(data = AllIMD, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllIMD, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="IMD")

####Stratified for unsuitable####
#Create new vars
Unsuitable<-merge(x=Unsuitable, y=Unique, by="patid", all.x=TRUE, all.y=TRUE)
Unsuitable[,2:15][is.na(Unsuitable[,2:15])]<-0

#By country
#Calculate totals

AllAny<-Unsuitable%>%
  group_by(country)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Unsuitable%>%
  group_by(country)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-country)

AllCountry<-cbind(AllActive, AllAny)

AllCountry$Year<-str_replace(AllCountry$ActYear, "Act", "")

AllCountry<-select(AllCountry, strat=country, Year, AnyCount, ActCount)
AllCountry$Group<-"Country"

CountryCI<-BinomCI(AllCountry$AnyCount, AllCountry$ActCount)*100

AllCountry<-cbind(AllCountry, CountryCI)

ggplot()+
  geom_line(data = AllCountry, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllCountry, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Country")

#By gender
AllAny<-Unsuitable%>%
  group_by(gender)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Unsuitable%>%
  group_by(gender)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-gender)

AllGender<-cbind(AllActive, AllAny)

AllGender$Year<-str_replace(AllGender$ActYear, "Act", "")

AllGender<-select(AllGender, strat=gender, Year, AnyCount, ActCount)
AllGender$Group<-"Sex"

GenderCI<-BinomCI(AllGender$AnyCount, AllGender$ActCount)*100

AllGender<-cbind(AllGender, GenderCI)

ggplot()+
  geom_line(data = AllGender, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllGender, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Sex")

#By ethnicity
AllAny<-Unsuitable%>%
  group_by(ethn_missing)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Unsuitable%>%
  group_by(ethn_missing)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-ethn_missing)

AllEthn<-cbind(AllActive, AllAny)

AllEthn$Year<-str_replace(AllEthn$ActYear, "Act", "")

AllEthn<-select(AllEthn, strat=ethn_missing, Year, AnyCount, ActCount)
AllEthn$Group<-"Ethnicity"

EthnCI<-BinomCI(AllEthn$AnyCount, AllEthn$ActCount)*100

AllEthn<-cbind(AllEthn, EthnCI)

ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllEthn, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Ethnicity")

#Without CI
ggplot()+
  geom_line(data = AllEthn, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Ethnicity")

####Now think about time varying ones####
####Age####

Age<-Unsuitable%>%
  mutate(Age0405 = case_when(2004-yob>=18 & 2004-yob<=30~ "18-30",
                             2004-yob>30 & 2004-yob<=40~ "31-40",
                             2004-yob>40 & 2004-yob<=50~ "41-50",
                             2004-yob>50 & 2004-yob<=60~ "51-60",
                             2004-yob>60 & 2004-yob<=70~ "61-70",
                             2004-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0506 = case_when(2005-yob>=18 & 2005-yob<=30~ "18-30",
                             2005-yob>30 & 2005-yob<=40~ "31-40",
                             2005-yob>40 & 2005-yob<=50~ "41-50",
                             2005-yob>50 & 2005-yob<=60~ "51-60",
                             2005-yob>60 & 2005-yob<=70~ "61-70",
                             2005-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0607 = case_when(2006-yob>=18 & 2006-yob<=30~ "18-30",
                             2006-yob>30 & 2006-yob<=40~ "31-40",
                             2006-yob>40 & 2006-yob<=50~ "41-50",
                             2006-yob>50 & 2006-yob<=60~ "51-60",
                             2006-yob>60 & 2006-yob<=70~ "61-70",
                             2006-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0708 = case_when(2007-yob>=18 & 2007-yob<=30~ "18-30",
                             2007-yob>30 & 2007-yob<=40~ "31-40",
                             2007-yob>40 & 2007-yob<=50~ "41-50",
                             2007-yob>50 & 2007-yob<=60~ "51-60",
                             2007-yob>60 & 2007-yob<=70~ "61-70",
                             2007-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0809 = case_when(2008-yob>=18 & 2008-yob<=30~ "18-30",
                             2008-yob>30 & 2008-yob<=40~ "31-40",
                             2008-yob>40 & 2008-yob<=50~ "41-50",
                             2008-yob>50 & 2008-yob<=60~ "51-60",
                             2008-yob>60 & 2008-yob<=70~ "61-70",
                             2008-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age0910 = case_when(2009-yob>=18 & 2009-yob<=30~ "18-30",
                             2009-yob>30 & 2009-yob<=40~ "31-40",
                             2009-yob>40 & 2009-yob<=50~ "41-50",
                             2009-yob>50 & 2009-yob<=60~ "51-60",
                             2009-yob>60 & 2009-yob<=70~ "61-70",
                             2009-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1011 = case_when(2010-yob>=18 & 2010-yob<=30~ "18-30",
                             2010-yob>30 & 2010-yob<=40~ "31-40",
                             2010-yob>40 & 2010-yob<=50~ "41-50",
                             2010-yob>50 & 2010-yob<=60~ "51-60",
                             2010-yob>60 & 2010-yob<=70~ "61-70",
                             2010-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1112 = case_when(2011-yob>=18 & 2011-yob<=30~ "18-30",
                             2011-yob>30 & 2011-yob<=40~ "31-40",
                             2011-yob>40 & 2011-yob<=50~ "41-50",
                             2011-yob>50 & 2011-yob<=60~ "51-60",
                             2011-yob>60 & 2011-yob<=70~ "61-70",
                             2011-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1213 = case_when(2012-yob>=18 & 2012-yob<=30~ "18-30",
                             2012-yob>30 & 2012-yob<=40~ "31-40",
                             2012-yob>40 & 2012-yob<=50~ "41-50",
                             2012-yob>50 & 2012-yob<=60~ "51-60",
                             2012-yob>60 & 2012-yob<=70~ "61-70",
                             2012-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1314 = case_when(2013-yob>=18 & 2013-yob<=30~ "18-30",
                             2013-yob>30 & 2013-yob<=40~ "31-40",
                             2013-yob>40 & 2013-yob<=50~ "41-50",
                             2013-yob>50 & 2013-yob<=60~ "51-60",
                             2013-yob>60 & 2013-yob<=70~ "61-70",
                             2013-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1415 = case_when(2014-yob>=18 & 2014-yob<=30~ "18-30",
                             2014-yob>30 & 2014-yob<=40~ "31-40",
                             2014-yob>40 & 2014-yob<=50~ "41-50",
                             2014-yob>50 & 2014-yob<=60~ "51-60",
                             2014-yob>60 & 2014-yob<=70~ "61-70",
                             2014-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1516 = case_when(2015-yob>=18 & 2015-yob<=30~ "18-30",
                             2015-yob>30 & 2015-yob<=40~ "31-40",
                             2015-yob>40 & 2015-yob<=50~ "41-50",
                             2015-yob>50 & 2015-yob<=60~ "51-60",
                             2015-yob>60 & 2015-yob<=70~ "61-70",
                             2015-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1617 = case_when(2016-yob>=18 & 2016-yob<=30~ "18-30",
                             2016-yob>30 & 2016-yob<=40~ "31-40",
                             2016-yob>40 & 2016-yob<=50~ "41-50",
                             2016-yob>50 & 2016-yob<=60~ "51-60",
                             2016-yob>60 & 2016-yob<=70~ "61-70",
                             2016-yob>70 ~ ">70",
                             TRUE ~ "HELP"),
         Age1718 = case_when(2017-yob>=18 & 2017-yob<=30~ "18-30",
                             2017-yob>30 & 2017-yob<=40~ "31-40",
                             2017-yob>40 & 2017-yob<=50~ "41-50",
                             2017-yob>50 & 2017-yob<=60~ "51-60",
                             2017-yob>60 & 2017-yob<=70~ "61-70",
                             2017-yob>70 ~ ">70",
                             TRUE ~ "HELP"))
         
#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAge<-subset(Age, Act0405==1)%>%
  select(patid, Age=Age0405, Any=Any0405)%>%
  group_by(Age)%>%
  summarise(ActCount=n(),AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
  NewFile<-subset(Age, Age[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Age=paste0("Age", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Age)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllAge<-rbind(AllAge, NewFile)
} 

AllAge$Age<-as.factor(AllAge$Age)
levels(AllAge$Age)<-c("18-30", "31-40", "41-50", "51-60", "61-70", ">70")

AllAge<-select(AllAge, strat=Age, Year, AnyCount, ActCount)
AllAge$Group<-"Age"

AllAgeCI<-BinomCI(AllAge$AnyCount, AllAge$ActCount)*100

AllAge<-cbind(AllAge, AllAgeCI)

ggplot()+
  geom_line(data = AllAge, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAge, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="Age")

####Diagnosis long####

Diag<-merge(x=Unsuitable, y=Diag3, by="patid", all.x=TRUE, All.y=FALSE)

DiagVar<-vars_select(names(Diag), starts_with('Diag', ignore.case = TRUE))

#Run it by year
AllDiag<-subset(Diag, Act0405==1)%>%
  select(patid, Diag=Diag0405, Any=Any0405)%>%
  group_by(Diag)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(DiagVar))) {
  NewFile<-subset(Diag, Diag[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Diag=paste0("Diag", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Diag)%>%
    summarise(ActCount=n(), AnyCount=sum(Any),Year=as.character(paste0(YearVar[i])))  
  AllDiag<-rbind(AllDiag, NewFile)
} 

AllDiag$Diag<-as.factor(AllDiag$Diag)

levels(AllDiag$Diag)<-c("Bipolar disorder", "Other psychoses", "Schizophrenia")
AllDiag$Diag <- factor(AllDiag$Diag, levels = c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

AllDiag<-select(AllDiag, strat=Diag, Year, AnyCount, ActCount)
AllDiag$Group<-"DiagTV"

AllDiagCI<-BinomCI(AllDiag$AnyCount, AllDiag$ActCount)*100

AllDiag<-cbind(AllDiag, AllDiagCI)

ggplot()+
  geom_line(data = AllDiag, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllDiag, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="SMI diagnosis")

####Antipsychotics####
Unsuitable<-merge(x=Unsuitable, y=APFields, by="patid", all.x=TRUE, all.y=FALSE)

AP<-Unsuitable%>%
  mutate(AP0405 = case_when(FirstAnyPres<="2005-03-31" ~ 1,
                            TRUE ~ 0),
         AP0506 = case_when(FirstAnyPres<="2006-03-31" ~ 1,
                            TRUE ~ 0),
         AP0607 = case_when(FirstAnyPres<="2007-03-31" ~ 1,
                            TRUE ~ 0),
         AP0708 = case_when(FirstAnyPres<="2008-03-31" ~ 1,
                            TRUE ~ 0), 
         AP0809 = case_when(FirstAnyPres<="2009-03-31"  ~ 1,
                            TRUE ~ 0),
         AP0910 = case_when(FirstAnyPres<="2010-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1011 = case_when(FirstAnyPres<="2011-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1112 = case_when(FirstAnyPres<="2012-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1213 = case_when(FirstAnyPres<="2013-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1314 = case_when(FirstAnyPres<="2014-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1415 = case_when(FirstAnyPres<="2015-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1516 = case_when(FirstAnyPres<="2016-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1617 = case_when(FirstAnyPres<="2017-03-31" ~ 1,
                            TRUE ~ 0), 
         AP1718 = case_when(FirstAnyPres<="2018-03-31" ~ 1,
                            TRUE ~ 0))

#Run it by year
YearVar<-(c("0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAP<-subset(AP, Act0405==1)%>%
  select(patid, AP=AP0405, Any=Any0405)%>%
  group_by(AP)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0405")

for (i in (2:length(YearVar))) {
  NewFile<-subset(AP, AP[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, AP=paste0("AP", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(AP)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllAP<-rbind(AllAP, NewFile)
} 

AllAP$AP<-as.factor(AllAP$AP)
levels(AllAP$AP)<-c("No", "Yes")

AllAP<-select(AllAP, strat=AP, Year, AnyCount, ActCount)
AllAP$Group<-"APTV"

AllAPCI<-BinomCI(AllAP$AnyCount, AllAP$ActCount)*100

AllAP<-cbind(AllAP, AllAPCI)

ggplot()+
  geom_line(data = AllAP, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllAP, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.3, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0,26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="On antipsychotics/mood stabilisers")

####For those with deprivation data####
AllAny<-Unsuitable%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-Unsuitable%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-FullIMD)

AllIMD<-cbind(AllActive, AllAny)

AllIMD$Year<-str_replace(AllIMD$ActYear, "Act", "")

AllIMD<-select(AllIMD, strat=FullIMD, Year, AnyCount, ActCount)
AllIMD$Group<-"IMD"

IMDCI<-BinomCI(AllIMD$AnyCount, AllIMD$ActCount)*100

AllIMD<-cbind(AllIMD, IMDCI)

ggplot()+
  geom_line(data = AllIMD, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllIMD, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 14), legend.title = element_text(colour = "black", size = 16))+
  scale_y_continuous(limits = c(0, 26), breaks = seq(0,26,by = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients exception reported", color="IMD")

####lOGISTIC REGRESSION#####
CPRD<-CPRD%>%
  mutate(YearAtEnd = case_when(end>="2004-04-01" & end<="2005-03-31" ~ "0405",
                               end>="2005-04-01" & end<="2006-03-31" ~ "0506",
                               end>="2006-04-01" & end<="2007-03-31" ~ "0607",
                               end>="2007-04-01" & end<="2008-03-31" ~ "0708",
                               end>="2008-04-01" & end<="2009-03-31" ~ "0809",
                               end>="2009-04-01" & end<="2010-03-31" ~ "0910",
                               end>="2010-04-01" & end<="2011-03-31" ~ "1011",
                               end>="2011-04-01" & end<="2012-03-31" ~ "1112",
                               end>="2012-04-01" & end<="2013-03-31" ~ "1213",
                               end>="2013-04-01" & end<="2014-03-31" ~ "1314",
                               end>="2014-04-01" & end<="2015-03-31" ~ "1415",
                               end>="2015-04-01" & end<="2016-03-31" ~ "1516",
                               end>="2016-04-01" & end<="2017-03-31" ~ "1617",
                               end>="2017-04-01" & end<="2018-03-31" ~ "1718",
                               TRUE ~ "HELP"))

table(CPRD$YearAtEnd)
CPRD$YearAtEnd<-as.factor(CPRD$YearAtEnd)
CPRD$TimeSinceReg<-CPRD$end-CPRD$crd

CPRD$LastDiag<-relevel(CPRD$LastDiag, ref = "bipolar")
CPRD$ethn_missing<-relevel(CPRD$ethn_missing, ref = "White")
CPRD$YearAtEnd<-relevel(CPRD$YearAtEnd, ref = "1718")

#All
RegAll<-glm(MH~AgeAtStart+gender+ethn_missing+country+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU+YearAtEnd, data=CPRD, family=binomial)

Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll<-RegAll%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAll.csv")

#Declined
RegDecline<-glm(EverDeclinedMH~AgeAtStart+gender+ethn_missing+country+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU+YearAtEnd, data=CPRD, family=binomial)

Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)

RegDecline<-RegDecline%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDecline.csv")
#Unsuitable
RegUnsuitable<-glm(EverUnsuitableMH~AgeAtStart+YearAtEnd+gender+ethn_missing+country+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU, data=CPRD, family=binomial)

Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable<-RegUnsuitable%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))

RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitable.csv")

###Unadj all###

#Age
RegAll<-glm(MH~AgeAtStart, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll<-RegAll%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAge.csv")

#gender
RegAll<-glm(MH~gender, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExReggender.csv")

#eth
RegAll<-glm(MH~ethn_missing, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegeth.csv")

#country
RegAll<-glm(MH~country, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegcountry.csv")

#LastDiag
RegAll<-glm(MH~LastDiag, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegLastDiag.csv")

#OtherQOF
RegAll<-glm(MH~OtherQOF, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegOtherQOF.csv")

#AnyPres
RegAll<-glm(MH~AnyPres, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAnyPres.csv")

#TimeSinceDiag
RegAll<-glm(MH~TimeSinceDiag, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegTimeSinceDiag.csv")

#TimeSinceReg
RegAll<-glm(MH~TimeSinceReg, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegTimeSinceReg.csv")

#FU
RegAll<-glm(MH~FU, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegFU.csv")

#YearAtEnd
RegAll<-glm(MH~YearAtEnd, data=CPRD, family=binomial)
Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegYearAtEnd.csv")

###Unadj declined###

#Age
RegDecline<-glm(EverDeclinedMH~AgeAtStart, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline<-RegDecline%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineAge.csv")

#gender
RegDecline<-glm(EverDeclinedMH~gender, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclinegender.csv")

#eth
RegDecline<-glm(EverDeclinedMH~ethn_missing, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineeth.csv")

#country
RegDecline<-glm(EverDeclinedMH~country, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclinecountry.csv")

#LastDiag
RegDecline<-glm(EverDeclinedMH~LastDiag, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineLastDiag.csv")

#OtherQOF
RegDecline<-glm(EverDeclinedMH~OtherQOF, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineOtherQOF.csv")

#AnyPres
RegDecline<-glm(EverDeclinedMH~AnyPres, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineAnyPres.csv")

#TimeSinceDiag
RegDecline<-glm(EverDeclinedMH~TimeSinceDiag, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineTimeSinceDiag.csv")

#TimeSinceReg
RegDecline<-glm(EverDeclinedMH~TimeSinceReg, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineTimeSinceReg.csv")

#FU
RegDecline<-glm(EverDeclinedMH~FU, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineFU.csv")

#YearAtEnd
RegDecline<-glm(EverDeclinedMH~YearAtEnd, data=CPRD, family=binomial)
Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineYearAtEnd.csv")

###unadj unsuitable###
#Age
RegUnsuitable<-glm(EverUnsuitableMH~AgeAtStart, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable<-RegUnsuitable%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableAge.csv")

#gender
RegUnsuitable<-glm(EverUnsuitableMH~gender, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitablegender.csv")

#eth
RegUnsuitable<-glm(EverUnsuitableMH~ethn_missing, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableeth.csv")

#country
RegUnsuitable<-glm(EverUnsuitableMH~country, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitablecountry.csv")

#LastDiag
RegUnsuitable<-glm(EverUnsuitableMH~LastDiag, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableLastDiag.csv")

#OtherQOF
RegUnsuitable<-glm(EverUnsuitableMH~OtherQOF, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableOtherQOF.csv")

#AnyPres
RegUnsuitable<-glm(EverUnsuitableMH~AnyPres, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableAnyPres.csv")

#TimeSinceDiag
RegUnsuitable<-glm(EverUnsuitableMH~TimeSinceDiag, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableTimeSinceDiag.csv")

#TimeSinceReg
RegUnsuitable<-glm(EverUnsuitableMH~TimeSinceReg, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableTimeSinceReg.csv")

#FU
RegUnsuitable<-glm(EverUnsuitableMH~FU, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableFU.csv")

#YearAtEnd
RegUnsuitable<-glm(EverUnsuitableMH~YearAtEnd, data=CPRD, family=binomial)
Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableYearAtEnd.csv")

####IMD Regression####
IMD<-subset(CPRD, !is.na(FullIMD))

#All
RegAll<-glm(MH~FullIMD+AgeAtStart+YearAtEnd+gender+ethn_missing+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU, data=IMD, family=binomial)

Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =IMD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)
RegAll<-RegAll%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAllIMD.csv")

#Declined
RegDecline<-glm(EverDeclinedMH~FullIMD+AgeAtStart+YearAtEnd+gender+ethn_missing+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU, data=IMD, family=binomial)

Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =IMD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)
RegDecline<-RegDecline%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))
RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDeclineIMD.csv")
#Unsuitable
RegUnsuitable<-glm(EverUnsuitableMH~FullIMD+AgeAtStart+YearAtEnd+gender+ethn_missing+LastDiag+OtherQOF+AnyPres+TimeSinceDiag+TimeSinceReg+FU, data=IMD, family=binomial)

Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =IMD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)
RegUnsuitable<-RegUnsuitable%>%
  mutate(estimate = case_when(term=="AgeAtStart" ~ estimate^10,
                              TRUE ~ estimate),
         `2.5 %` = case_when(term=="AgeAtStart" ~ `2.5 %`^10,
                             TRUE ~ `2.5 %`),
         `97.5 %` = case_when(term=="AgeAtStart" ~ `97.5 %`^10,
                              TRUE ~ `97.5 %`))
RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitableIMD.csv")

####By year####
#2004-2011
which(colnames(CPRD) == "Act0405")
which(colnames(CPRD) == "Act1011")
which(colnames(CPRD) == "Act1314")
which(colnames(CPRD) == "Act1718")

CPRD$Start04<-pmax(CPRD$start, as.Date("2004-04-01"))
CPRD$End04<-pmin(CPRD$end, as.Date("2011-03-31"))

CPRD$Start11<-pmax(CPRD$start, as.Date("2011-04-01"))
CPRD$End11<-pmin(CPRD$end, as.Date("2014-03-31"))

CPRD$Start14<-pmax(CPRD$start, as.Date("2014-04-01"))
CPRD$End14<-pmin(CPRD$end, as.Date("2018-03-31"))

CPRDMod<-CPRD%>%
  ungroup()%>%
  rowwise()%>%
  mutate(Act04 = sum(c_across(129:135)),
         Act11 = sum(c_across(136:138)),
         Act18 = sum(c_across(139:142)))

length(which(CPRDMod$Act04>=1 | CPRDMod$Act11>=1|CPRDMod$Act18>=1))

CPRD2004<-subset(CPRDMod, Act04>=1)
CPRD2011<-subset(CPRDMod, Act11>=1)
CPRD2014<-subset(CPRDMod, Act18>=1)

#year of end of follow up
CPRD2004<-CPRD2004%>%
  mutate(Year04 = case_when(End04>="2004-04-01" & End04<="2005-03-31" ~ "0405",
                            End04>="2005-04-01" & End04<="2006-03-31" ~ "0506",
                            End04>="2006-04-01" & End04<="2007-03-31" ~ "0607",
                            End04>="2007-04-01" & End04<="2008-03-31" ~ "0708",
                            End04>="2008-04-01" & End04<="2009-03-31" ~ "0809",
                            End04>="2009-04-01" & End04<="2010-03-31" ~ "0910",
                            End04>="2010-04-01" & End04<="2011-03-31" ~ "1011",
                            TRUE ~ "HELP"))
table(CPRD2004$Year04)

CPRD2011<-CPRD2011%>%
  mutate(Year11 = case_when(End11>="2011-04-01" & End11<="2012-03-31" ~ "1112",
                            End11>="2012-04-01" & End11<="2013-03-31" ~ "1213",
                            End11>="2013-04-01" & End11<="2014-03-31" ~ "1314",
                            TRUE ~ "HELP"))
table(CPRD2011$Year11)

CPRD2014<-CPRD2014%>%
  mutate(Year14 = case_when(End14>="2014-04-01" & End14<="2015-03-31" ~ "1415",
                            End14>="2015-04-01" & End14<="2016-03-31" ~ "1516",
                            End14>="2016-04-01" & End14<="2017-03-31" ~ "1617",
                            End14>="2017-04-01" & End14<="2018-03-31" ~ "1718",
                            TRUE ~ "HELP"))
table(CPRD2014$Year14)

CPRD2004$Year04<-as.factor(CPRD2004$Year04)
CPRD2011$Year11<-as.factor(CPRD2011$Year11)
CPRD2014$Year14<-as.factor(CPRD2014$Year14)

CPRD2004$Year04 <- relevel(CPRD2004$Year04, ref = "1011")
CPRD2011$Year11 <- relevel(CPRD2011$Year11, ref = "1314")
CPRD2014$Year14 <- relevel(CPRD2014$Year14, ref = "1718")

CPRD$AgeAtStart04<-as.numeric(year(CPRD$Start04)-CPRD$yob)
CPRD$AgeAtStart11<-as.numeric(year(CPRD$Start11)-CPRD$yob)
CPRD$AgeAtStart14<-as.numeric(year(CPRD$Start14)-CPRD$yob)

CPRD$TimeSinceReg04<-CPRD$End04-CPRD$crd
CPRD$TimeSinceReg11<-CPRD$End11-CPRD$crd
CPRD$TimeSinceReg14<-CPRD$End14-CPRD$crd

CPRD$TimeSinceDiag04<-CPRD$End04-CPRD$FirstSMIDate
CPRD$TimeSinceDiag11<-CPRD$End11-CPRD$FirstSMIDate
CPRD$TimeSinceDiag14<-CPRD$End14-CPRD$FirstSMIDate

CPRD$TimeSinceReg04<-CPRD$End04-CPRD$crd
CPRD$TimeSinceReg11<-CPRD$End11-CPRD$crd
CPRD$TimeSinceReg14<-CPRD$End14-CPRD$crd

CPRD$FU04<-CPRD$End04-CPRD$Start04
CPRD$FU11<-CPRD$End11-CPRD$Start11
CPRD$FU14<-CPRD$End14-CPRD$Start14

#In period variables
#Exceptions to include non-QOF "physical health check exemption"
load("VariableExtracts/CPRD2018/MergedObs/SMIEx.Rdata")
PatExAll<-subset(PatExAll, patid %in% CPRD$patid & Group=="MH")

Dates<-select(CPRD, patid, start, end, Start04, End04, End11, Start11, End14, Start14, yob)
PatExAll<-merge(x=PatExAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatExAll<-subset(PatExAll, eventdate<=end & eventdate>yob)

Exempt2004<-subset(PatExAll, eventdate<=End04)
Exempt2011<-subset(PatExAll, eventdate<=End11)
Exempt2014<-subset(PatExAll, eventdate<=End14)

CPRD$EverMHEx2004<-0
CPRD$EverMHEx2004[CPRD$patid %in% Exempt2004$patid]<-1

CPRD$EverMHEx2011<-0
CPRD$EverMHEx2011[CPRD$patid %in% Exempt2011$patid]<-1

CPRD$EverMHEx2014<-0
CPRD$EverMHEx2014[CPRD$patid %in% Exempt2014$patid]<-1

#Other QOF
load("VariableExtracts/CPRD2018/MergedObs/OtherQOF.Rdata")
PatOtherQOFAll<-subset(PatOtherQOFAll, patid %in% CPRD$patid)

PatOtherQOFAll<-merge(x=PatOtherQOFAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatOtherQOFAll<-subset(PatOtherQOFAll, eventdate<=end & eventdate>yob)

OtherQOF2004<-subset(PatOtherQOFAll, eventdate<=End04)
OtherQOF2011<-subset(PatOtherQOFAll, eventdate<=End11)
OtherQOF2014<-subset(PatOtherQOFAll, eventdate<=End14)

CPRD$EverOtherQOF2004<-0
CPRD$EverOtherQOF2004[CPRD$patid %in% OtherQOF2004$patid]<-1

CPRD$EverOtherQOF2011<-0
CPRD$EverOtherQOF2011[CPRD$patid %in% OtherQOF2011$patid]<-1

CPRD$EverOtherQOF2014<-0
CPRD$EverOtherQOF2014[CPRD$patid %in% OtherQOF2014$patid]<-1

#Antipsych
load("VariableExtracts/CPRD2018/MergedObs/APProd.Rdata")
PatAPAll<-subset(PatAPAll, patid %in% CPRD$patid)
PatAPAll<-merge(x=PatAPAll, y=Dates, by="patid", all.x=TRUE, all.y=FALSE)
PatAPAll<-subset(PatAPAll, eventdate>=start & eventdate<=end)

AP2004<-subset(PatAPAll, eventdate<=End04)
AP2011<-subset(PatAPAll, eventdate<=End11)
AP2014<-subset(PatAPAll, eventdate<=End14)

CPRD$AP2004<-0
CPRD$AP2004[CPRD$patid %in% AP2004$patid]<-1

CPRD$AP2011<-0
CPRD$AP2011[CPRD$patid %in% AP2011$patid]<-1

CPRD$AP2014<-0
CPRD$AP2014[CPRD$patid %in% AP2014$patid]<-1


##2004
#All
RegAll<-glm(MH~AgeAtStart04+End04+gender+ethn_missing+country+LastDiag+EverOtherQOF2004+AP2004+TimeSinceDiag04+TimeSinceReg04+FU04, data=CPRD, family=binomial)

Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAll04.csv")

#Declined
RegDecline<-glm(EverDeclinedMH~AgeAtStart04+End04+gender+ethn_missing+country+LastDiag+EverOtherQOF2004+AP2004+TimeSinceDiag04+TimeSinceReg04+FU04, data=CPRD, family=binomial)

Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)

RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDecline04.csv")

#Unsuitable
RegUnsuitable<-glm(EverUnsuitableMH~AgeAtStart04+End04+gender+ethn_missing+country+LastDiag+EverOtherQOF2004+AP2004+TimeSinceDiag04+TimeSinceReg04+FU04, data=CPRD, family=binomial)

Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)

RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitable04.csv")

##2011
#All
RegAll<-glm(MH~AgeAtStart11+End11+gender+ethn_missing+country+LastDiag+EverOtherQOF2011+AP2011+TimeSinceDiag11+TimeSinceReg11+FU11, data=CPRD, family=binomial)

Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAll11.csv")

#Declined
RegDecline<-glm(EverDeclinedMH~AgeAtStart11+End11+gender+ethn_missing+country+LastDiag+EverOtherQOF2011+AP2011+TimeSinceDiag11+TimeSinceReg11+FU11, data=CPRD, family=binomial)

Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)

RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDecline11.csv")

#Unsuitable
RegUnsuitable<-glm(EverUnsuitableMH~AgeAtStart11+End11+gender+ethn_missing+country+LastDiag+EverOtherQOF2011+AP2011+TimeSinceDiag11+TimeSinceReg11+FU11, data=CPRD, family=binomial)

Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)

RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitable11.csv")

##2014
#All
RegAll<-glm(MH~AgeAtStart14+End14+gender+ethn_missing+country+LastDiag+EverOtherQOF2014+AP2014+TimeSinceDiag14+TimeSinceReg14+FU14, data=CPRD, family=binomial)

Results<-tidy(RegAll, exponentiate=TRUE)
Sand<-vcovCL(RegAll, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegAll, vcov=Sand)))
RegAll<-cbind(Results, CIFin)

RegAll$Result<-paste0(format(round(RegAll$estimate, 2), nsmall=2), " (", format(round(RegAll$`2.5 %`, 2), nsmall=2), "-", format(round(RegAll$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegAll, file="DatamindCPRD/Outputs/ExRegAll14.csv")

#Declined
RegDecline<-glm(EverDeclinedMH~AgeAtStart14+End14+gender+ethn_missing+country+LastDiag+EverOtherQOF2014+AP2014+TimeSinceDiag14+TimeSinceReg14+FU14, data=CPRD, family=binomial)

Results<-tidy(RegDecline, exponentiate=TRUE)
Sand<-vcovCL(RegDecline, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegDecline, vcov=Sand)))
RegDecline<-cbind(Results, CIFin)

RegDecline$Result<-paste0(format(round(RegDecline$estimate, 2), nsmall=2), " (", format(round(RegDecline$`2.5 %`, 2), nsmall=2), "-", format(round(RegDecline$`97.5 %`, 2), nsmall=2), ")")

write.csv(RegDecline, file="DatamindCPRD/Outputs/ExRegDecline14.csv")

#Unsuitable
RegUnsuitable<-glm(EverUnsuitableMH~AgeAtStart14+End14+gender+ethn_missing+country+LastDiag+EverOtherQOF2014+AP2014+TimeSinceDiag14+TimeSinceReg14+FU14, data=CPRD, family=binomial)

Results<-tidy(RegUnsuitable, exponentiate=TRUE)
Sand<-vcovCL(RegUnsuitable, cluster =CPRD$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(RegUnsuitable, vcov=Sand)))
RegUnsuitable<-cbind(Results, CIFin)

RegUnsuitable$Result<-paste0(format(round(RegUnsuitable$estimate, 2), nsmall=2), " (", format(round(RegUnsuitable$`2.5 %`, 2), nsmall=2), "-", format(round(RegUnsuitable$`97.5 %`, 2), nsmall=2), ")")
write.csv(RegUnsuitable, file="DatamindCPRD/Outputs/ExRegUnsuitable14.csv")



