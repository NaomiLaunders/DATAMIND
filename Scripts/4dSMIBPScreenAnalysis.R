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

####Load SMI cohort####

load("DatamindCPRD/Data/MainAnalysis.Rdata")

####Set active CPRD cohort####
CPRDActive<-CPRD%>%
  select(patid, Act0001, Act0102, Act0203, Act0304, Act0405, Act0506, Act0607, Act0708, Act0809, Act0910, Act1011, Act1112, Act1213, Act1314, Act1415, Act1516, Act1617, Act1718)
CPRDActive[,c(2:19)]<-lapply(CPRDActive[,c(2:19)], as.character)
CPRDActive[,c(2:19)]<-lapply(CPRDActive[,c(2:19)], as.numeric)

cbPalette <- c("#999999", "#CC79A7", "#E69F00","#009E73", "#F0E442", "#0072B2")

####BP####
load("DatamindCPRD/Data/BPAnalysis.Rdata")

Any<-BP%>%
  mutate(Year = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                          eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                          eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                          eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                          eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
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

AnyBP<-Any%>%
  group_by(patid, Year)%>%
  summarise(Count=n())%>%
  arrange(Year)%>%
  pivot_wider(names_from=Year, names_prefix = "Any",  values_from = Count, values_fill = 0)

#Merge to active cohort

All<-merge(x=CPRDActive, y=AnyBP, by="patid", all.x=TRUE, all.y=TRUE)

#Change to binary and set to 0 if not active
ActVar<-vars_select(names(CPRDActive), starts_with('Act', ignore.case = TRUE))
BPVar<-vars_select(names(AnyBP), starts_with('Any', ignore.case = TRUE))

for (i in (1:length(BPVar))) {
  All[paste0(BPVar[i])]<-if_else(All[paste0(ActVar[i])]==0, 0, if_else(All[paste0(ActVar[i])]==1 & All[paste0(BPVar[i])]>0, 1, 0))
  All[paste0(BPVar[i])][is.na(All[paste0(BPVar[i])])]<-0
}

####Stratification####

#Create new vars
NewVars<-select(CPRD, patid, country, LastDiag, gender, yob, FirstSpecificQOF, pracid, ethn_missing, EverMHEx, EverCholEx, EverQOFExempt, FullIMD)
All<-merge(x=All, y=NewVars, by="patid", all.x=TRUE, all.y=TRUE)

#Just double check
length(which(All$Any1718==1 & All$Act1718==0))

#By country
#Calculate totals

AllAny<-All%>%
  group_by(country)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
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

AllCountry$Year<-case_when(AllCountry$Year=="0001" ~ "2000-2001",
                           AllCountry$Year=="0102" ~ "2001-2002",
                           AllCountry$Year=="0203" ~ "2002-2003",
                           AllCountry$Year=="0304" ~ "2003-2004",
                           AllCountry$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Country")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)
  

ggsave("DatamindCPRD/Outputs/CountryBP.PNG", width = 10,   height = 10)

#By gender
AllAny<-All%>%
  group_by(gender)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(gender)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-gender)

AllSex<-cbind(AllActive, AllAny)

AllSex$Year<-str_replace(AllSex$ActYear, "Act", "")

AllSex<-select(AllSex, strat=gender, Year, AnyCount, ActCount)
AllSex$Group<-"Sex"

SexCI<-BinomCI(AllSex$AnyCount, AllSex$ActCount)*100

AllSex<-cbind(AllSex, SexCI)

AllSex$Year<-case_when(AllSex$Year=="0001" ~ "2000-2001",
                       AllSex$Year=="0102" ~ "2001-2002",
                       AllSex$Year=="0203" ~ "2002-2003",
                       AllSex$Year=="0304" ~ "2003-2004",
                       AllSex$Year=="0405" ~ "2004-2005",
                       AllSex$Year=="0506" ~ "2005-2006",
                       AllSex$Year=="0607" ~ "2006-2007",
                       AllSex$Year=="0708" ~ "2007-2008",
                       AllSex$Year=="0809" ~ "2008-2009",
                       AllSex$Year=="0910" ~ "2009-2010",
                       AllSex$Year=="1011" ~ "2010-2011",
                       AllSex$Year=="1112" ~ "2011-2012",
                       AllSex$Year=="1213" ~ "2012-2013",
                       AllSex$Year=="1314" ~ "2013-2014",
                       AllSex$Year=="1415" ~ "2014-2015",
                       AllSex$Year=="1516" ~ "2015-2016",
                       AllSex$Year=="1617" ~ "2016-2017",
                       AllSex$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllSex, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllSex, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Sex")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/GenderBP.PNG", width = 10,   height = 10)

#By ethnicity
AllAny<-All%>%
  group_by(ethn_missing)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(ethn_missing)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-ethn_missing)

AllEthn<-cbind(AllActive, AllAny)

AllEthn$Year<-str_replace(AllEthn$ActYear, "Act", "")

AllEthn$Ethnicity<-factor(AllEthn$ethn_missing, levels = c("Asian", "Black", "Mixed", "White", "Other", "Missing"))

AllEthn<-select(AllEthn, strat=Ethnicity, Year, AnyCount, ActCount)
AllEthn$Group<-"Ethnicity"

EthnCI<-BinomCI(AllEthn$AnyCount, AllEthn$ActCount)*100

AllEthn<-cbind(AllEthn, EthnCI)

AllEthn$Year<-case_when(AllEthn$Year=="0001" ~ "2000-2001",
                        AllEthn$Year=="0102" ~ "2001-2002",
                        AllEthn$Year=="0203" ~ "2002-2003",
                        AllEthn$Year=="0304" ~ "2003-2004",
                        AllEthn$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Ethnicity")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/EthnicityBP.PNG", width = 10,   height = 10)

#By practice
AllAny<-All%>%
  group_by(pracid)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(pracid)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-pracid)

AllPrac<-cbind(AllActive, AllAny)

AllPrac$Anyprop<-AllPrac$AnyCount/AllPrac$ActCount*100
AllPrac$Year<-str_replace(AllPrac$ActYear, "Act", "")

#Limit to practices with at least 10 active patients
PracIDToDrop<-subset(AllPrac, ActCount<10)
AllPrac<-subset(AllPrac, !pracid %in% PracIDToDrop$pracid)
length(unique(AllPrac$pracid))
length(unique(CPRD$pracid))

AllPrac$Year<-case_when(AllPrac$Year=="0001" ~ "2000-2001",
                        AllPrac$Year=="0102" ~ "2001-2002",
                        AllPrac$Year=="0203" ~ "2002-2003",
                        AllPrac$Year=="0304" ~ "2003-2004",
                        AllPrac$Year=="0405" ~ "2004-2005",
                        AllPrac$Year=="0506" ~ "2005-2006",
                        AllPrac$Year=="0607" ~ "2006-2007",
                        AllPrac$Year=="0708" ~ "2007-2008",
                        AllPrac$Year=="0809" ~ "2008-2009",
                        AllPrac$Year=="0910" ~ "2009-2010",
                        AllPrac$Year=="1011" ~ "2010-2011",
                        AllPrac$Year=="1112" ~ "2011-2012",
                        AllPrac$Year=="1213" ~ "2012-2013",
                        AllPrac$Year=="1314" ~ "2013-2014",
                        AllPrac$Year=="1415" ~ "2014-2015",
                        AllPrac$Year=="1516" ~ "2015-2016",
                        AllPrac$Year=="1617" ~ "2016-2017",
                        AllPrac$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllPrac, stat = "identity", aes(x=Year, y=Anyprop, color=pracid, group=pracid))+
  theme_classic()+
  guides(color="none")+
  guides(fill="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/PracticeBP.PNG", width = 10,   height = 10)

#By diag
AllAny<-All%>%
  group_by(LastDiag)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(LastDiag)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-LastDiag)

AllDiag1<-cbind(AllActive, AllAny)

AllDiag1$Year<-str_replace(AllDiag1$ActYear, "Act", "")

AllDiag1$LastDiag<-as.factor(AllDiag1$LastDiag)
levels(AllDiag1$LastDiag)<-c("Bipolar disorder", "Other psychoses", "Schizophrenia")
AllDiag1$LastDiag <- factor(AllDiag1$LastDiag, levels = c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

AllDiag1<-select(AllDiag1, strat=LastDiag, Year, AnyCount, ActCount)
AllDiag1$Group<-"diag"

DiagCI<-BinomCI(AllDiag1$AnyCount, AllDiag1$ActCount)*100

AllDiag1<-cbind(AllDiag1, DiagCI)

ggplot()+
  geom_line(data = AllDiag1, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllDiag1, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients receiving BP screening", color="SMI diagnosis")

#Ever exempt
All$MHBin<-as.factor(All$EverMHEx)

AllAny<-All%>%
  group_by(MHBin)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(MHBin)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-MHBin)

AllEx1<-cbind(AllActive, AllAny)

AllEx1$Year<-str_replace(AllEx1$ActYear, "Act", "")

AllEx1<-select(AllEx1, strat=MHBin, Year, AnyCount, ActCount)
AllEx1$Group<-"Exempt"

ExCI<-BinomCI(AllEx1$AnyCount, AllEx1$ActCount)*100

AllEx1<-cbind(AllEx1, ExCI)

ggplot()+
  geom_line(data = AllEx1, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllEx1, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients receiving BP screening", color="Ever exempt")

####Other QOF ever####
All$OtherQOF<-0
All$OtherQOF[!is.na(All$FirstSpecificQOF)]<-1
All$OtherQOF<-as.factor(All$OtherQOF)

AllAny<-All%>%
  group_by(OtherQOF)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
  group_by(OtherQOF)%>%
  select(starts_with("Act"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Act"), names_to = "ActYear", values_to = "ActCount")%>%
  select(-OtherQOF)

AllQOF1<-cbind(AllActive, AllAny)

AllQOF1$Year<-str_replace(AllQOF1$ActYear, "Act", "")

AllQOF1<-select(AllQOF1, strat=OtherQOF, Year, AnyCount, ActCount)
AllQOF1$Group<-"QOF"

QOFCI<-BinomCI(AllQOF1$AnyCount, AllQOF1$ActCount)*100

AllQOF1<-cbind(AllQOF1, QOFCI)

ggplot()+
  geom_line(data = AllQOF1, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllQOF1, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(fill="none")+theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  labs( x="Year", y = "Percentage of patients receiving BP screening", color="Ever on another QOF register")

####Now think about time varying ones####
####Age####

Age<-All%>%
  mutate(Age0001 = case_when(2000-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0102 = case_when(2001-yob>=40~ 1,
                             TRUE ~ 0),
         Age0203 = case_when(2002-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0304 = case_when(2003-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0405 = case_when(2004-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0506 = case_when(2005-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0607 = case_when(2006-yob>=40 ~ 1,
                             TRUE ~ 0),
         Age0708 = case_when(2007-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age0809 = case_when(2008-yob>=40  ~ 1,
                             TRUE ~ 0), 
         Age0910 = case_when(2009-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1011 = case_when(2010-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1112 = case_when(2011-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1213 = case_when(2012-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1314 = case_when(2013-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1415 = case_when(2014-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1516 = case_when(2015-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1617 = case_when(2016-yob>=40 ~ 1,
                             TRUE ~ 0), 
         Age1718 = case_when(2017-yob>=40 ~ 1,
                             TRUE ~ 0))

#Run it by year


YearVar<-(c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAge<-subset(Age, Act0001==1)%>%
  select(patid, Age=Age0001, Any=Any0001)%>%
  group_by(Age)%>%
  summarise(ActCount=n(), AnyCount=sum(Any),Year="0001")

for (i in (2:length(YearVar))) {
  NewFile<-subset(Age, Age[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Age=paste0("Age", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Age)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllAge<-rbind(AllAge, NewFile)
} 

AllAge$Age<-as.factor(AllAge$Age)
levels(AllAge$Age)<-c("Under 40", "40+")

AllAge<-select(AllAge, strat=Age, Year, AnyCount, ActCount)
AllAge$Group<-"Age"

AllAgeCI<-BinomCI(AllAge$AnyCount, AllAge$ActCount)*100

AllAge<-cbind(AllAge, AllAgeCI)

AllAge$Year<-case_when(AllAge$Year=="0001" ~ "2000-2001",
                       AllAge$Year=="0102" ~ "2001-2002",
                       AllAge$Year=="0203" ~ "2002-2003",
                       AllAge$Year=="0304" ~ "2003-2004",
                       AllAge$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Age")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/AgeBP.PNG", width = 10,   height = 10)

####Exempt####
load("VariableExtracts/CPRD2018/MergedObs/SMIEx.Rdata")

PatExAll<-subset(PatExAll, patid %in% CPRD$patid)

length(unique(PatExAll$patid))

Exempt<-subset(PatExAll, eventdate>="2000-04-01" & eventdate<="2018-03-31")
length(unique(Exempt$patid))

table(Exempt$Group)
table(Exempt$Type)
Exempt$Type[Exempt$Type=="MH"]<-"Unsuitable"

Exempt<-Exempt%>%
  subset(Group=="MH")%>%
  group_by(patid)%>%
  mutate(ExYear = case_when(eventdate>="2000-04-01" & eventdate<="2001-03-31" ~ "0001",
                            eventdate>="2001-04-01" & eventdate<="2002-03-31" ~ "0102",
                            eventdate>="2002-04-01" & eventdate<="2003-03-31" ~ "0203",
                            eventdate>="2003-04-01" & eventdate<="2004-03-31" ~ "0304",
                            eventdate>="2004-04-01" & eventdate<="2005-03-31" ~ "0405",
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

Exempt<-Exempt%>%
  group_by(patid, ExYear)%>%
  summarise(Count=n())%>%
  arrange(ExYear)%>%
  pivot_wider(names_from=ExYear, names_prefix = "Ex",  values_from = Count, values_fill = 0)

Exempt<-merge(x=All, y=Exempt, by="patid", all.x=TRUE, All.y=FALSE)

ExVar<-vars_select(names(Exempt), starts_with('Ex', ignore.case = TRUE))

for (i in (1:length(ExVar))) {
  Exempt[paste0(ExVar[i])]<-if_else(Exempt[paste0(ExVar[i])]==0, 0, if_else(Exempt[paste0(ExVar[i])]>0 & Exempt[paste0(ActVar[i])]==1, 1, 0))
  Exempt[paste0(ExVar[i])][is.na(Exempt[paste0(ExVar[i])])]<-0
}

#Run it by year
AllExempt<-subset(Exempt, Act0001==1)%>%
  select(patid, Exempt=Ex0001, Any=Any0001)%>%
  group_by(Exempt)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0001")

for (i in (2:length(YearVar))) {
  NewFile<-subset(Exempt, Exempt[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Exempt=paste0("Ex", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Exempt)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllExempt<-rbind(AllExempt, NewFile)
} 

AllExempt$Exempt<-as.factor(AllExempt$Exempt)
levels(AllExempt$Exempt)<-c("Not Exempt", "Exempt")

AllExempt<-select(AllExempt, strat=Exempt, Year, AnyCount, ActCount)
AllExempt$Group<-"ExemptTV"

AllExemptCI<-BinomCI(AllExempt$AnyCount, AllExempt$ActCount)*100

AllExempt<-cbind(AllExempt, AllExemptCI)

AllExempt$Year<-case_when(AllExempt$Year=="0001" ~ "2000-2001",
                          AllExempt$Year=="0102" ~ "2001-2002",
                          AllExempt$Year=="0203" ~ "2002-2003",
                          AllExempt$Year=="0304" ~ "2003-2004",
                          AllExempt$Year=="0405" ~ "2004-2005",
                          AllExempt$Year=="0506" ~ "2005-2006",
                          AllExempt$Year=="0607" ~ "2006-2007",
                          AllExempt$Year=="0708" ~ "2007-2008",
                          AllExempt$Year=="0809" ~ "2008-2009",
                          AllExempt$Year=="0910" ~ "2009-2010",
                          AllExempt$Year=="1011" ~ "2010-2011",
                          AllExempt$Year=="1112" ~ "2011-2012",
                          AllExempt$Year=="1213" ~ "2012-2013",
                          AllExempt$Year=="1314" ~ "2013-2014",
                          AllExempt$Year=="1415" ~ "2014-2015",
                          AllExempt$Year=="1516" ~ "2015-2016",
                          AllExempt$Year=="1617" ~ "2016-2017",
                          AllExempt$Year=="1718" ~ "2017-2018")

AllExempt<-subset(AllExempt, Year!="2000-2001" & Year!= "2001-2002" & Year!="2002-2003" & Year!="2003-2004")

ggplot()+
  geom_line(data = AllExempt, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllExempt, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Exception reported")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=7.5, linetype='dotted')+
  annotate("text", x = 7.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/ExceptionsBP.PNG", width = 10,   height = 10)

####QOF####

Other<-All%>%
  mutate(Other0001 = case_when(FirstSpecificQOF<="2001-03-31" ~ 1,
                               TRUE ~ 0),
         Other0102 = case_when(FirstSpecificQOF<="2002-03-31"~ 1,
                               TRUE ~ 0),
         Other0203 = case_when(FirstSpecificQOF<="2003-03-31" ~ 1,
                               TRUE ~ 0),
         Other0304 = case_when(FirstSpecificQOF<="2004-03-31" ~ 1,
                               TRUE ~ 0),
         Other0405 = case_when(FirstSpecificQOF<="2005-03-31" ~ 1,
                               TRUE ~ 0),
         Other0506 = case_when(FirstSpecificQOF<="2006-03-31" ~ 1,
                               TRUE ~ 0),
         Other0607 = case_when(FirstSpecificQOF<="2007-03-31" ~ 1,
                               TRUE ~ 0),
         Other0708 = case_when(FirstSpecificQOF<="2008-03-31" ~ 1,
                               TRUE ~ 0), 
         Other0809 = case_when(FirstSpecificQOF<="2009-03-31"  ~ 1,
                               TRUE ~ 0),
         Other0910 = case_when(FirstSpecificQOF<="2010-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1011 = case_when(FirstSpecificQOF<="2011-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1112 = case_when(FirstSpecificQOF<="2012-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1213 = case_when(FirstSpecificQOF<="2013-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1314 = case_when(FirstSpecificQOF<="2014-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1415 = case_when(FirstSpecificQOF<="2015-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1516 = case_when(FirstSpecificQOF<="2016-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1617 = case_when(FirstSpecificQOF<="2017-03-31" ~ 1,
                               TRUE ~ 0), 
         Other1718 = case_when(FirstSpecificQOF<="2018-03-31" ~ 1,
                               TRUE ~ 0))

#Run it by year
YearVar<-(c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllQOF<-subset(Other, Act0001==1)%>%
  select(patid, Other=Other0001, Any=Any0001)%>%
  group_by(Other)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0001")

for (i in (2:length(YearVar))) {
  NewFile<-subset(Other, Other[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Other=paste0("Other", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Other)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllQOF<-rbind(AllQOF, NewFile)
} 

AllQOF$OtherQOF<-as.factor(AllQOF$Other)
levels(AllQOF$OtherQOF)<-c("No", "Yes")

AllQOF<-select(AllQOF, strat=OtherQOF, Year, AnyCount, ActCount)
AllQOF$Group<-"QOFTV"

AllQOFCI<-BinomCI(AllQOF$AnyCount, AllQOF$ActCount)*100

AllQOF<-cbind(AllQOF, AllQOFCI)

AllQOF$Year<-case_when(AllQOF$Year=="0001" ~ "2000-2001",
                       AllQOF$Year=="0102" ~ "2001-2002",
                       AllQOF$Year=="0203" ~ "2002-2003",
                       AllQOF$Year=="0304" ~ "2003-2004",
                       AllQOF$Year=="0405" ~ "2004-2005",
                       AllQOF$Year=="0506" ~ "2005-2006",
                       AllQOF$Year=="0607" ~ "2006-2007",
                       AllQOF$Year=="0708" ~ "2007-2008",
                       AllQOF$Year=="0809" ~ "2008-2009",
                       AllQOF$Year=="0910" ~ "2009-2010",
                       AllQOF$Year=="1011" ~ "2010-2011",
                       AllQOF$Year=="1112" ~ "2011-2012",
                       AllQOF$Year=="1213" ~ "2012-2013",
                       AllQOF$Year=="1314" ~ "2013-2014",
                       AllQOF$Year=="1415" ~ "2014-2015",
                       AllQOF$Year=="1516" ~ "2015-2016",
                       AllQOF$Year=="1617" ~ "2016-2017",
                       AllQOF$Year=="1718" ~ "2017-2018")

ggplot()+
  geom_line(data = AllQOF, stat = "identity", aes(x=Year, y=est, color=strat, group=strat))+
  geom_ribbon(data = AllQOF, aes(x=Year, ymin = lwr.ci, ymax = upr.ci, group=strat, fill=strat), alpha=0.5)+
  theme_classic()+
  guides(color="none")+
  theme(legend.position = c(0.2, 0.93),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="On another QOF register")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/QOFBP.PNG", width = 10,   height = 10)


####Diagnosis long####
load("VariableExtracts/CPRD2018/MergedObs/SMIObs.Rdata")
SMIObs<-subset(SMIObs, patid %in% CPRD$patid)

Diag<-SMIObs%>%
  group_by(patid)%>%
  mutate(DiagYear = case_when(SMIDate<"2000-04-01" ~ "0000",
                              SMIDate>="2000-04-01" & SMIDate<="2001-03-31" ~ "0001",
                              SMIDate>="2001-04-01" & SMIDate<="2002-03-31" ~ "0102",
                              SMIDate>="2002-04-01" & SMIDate<="2003-03-31" ~ "0203",
                              SMIDate>="2003-04-01" & SMIDate<="2004-03-31" ~ "0304",
                              SMIDate>="2004-04-01" & SMIDate<="2005-03-31" ~ "0405",
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
                              TRUE ~ NA_character_))
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
  pivot_longer(cols= c(2:20), names_to = "DiagYear", values_to = "diag")%>%
  group_by(patid)%>%
  mutate(diag=case_when(is.na(diag) ~ lag(diag),
                        TRUE ~ diag))%>%
  fill(diag, .direction=c("down"))

Diag3<-Diag2%>%
  subset(DiagYear!="Diag0000")%>%
  pivot_wider(names_from=DiagYear, values_from = diag)

Diag<-merge(x=All, y=Diag3, by="patid", all.x=TRUE, All.y=FALSE)

DiagVar<-vars_select(names(Diag), starts_with('Diag', ignore.case = TRUE))

#Run it by year
AllDiag<-subset(Diag, Act0001==1)%>%
  select(patid, Diag=Diag0001, Any=Any0001)%>%
  group_by(Diag)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0001")

for (i in (2:length(DiagVar))) {
  NewFile<-subset(Diag, Diag[paste0("Act", YearVar[i])]==1)
  NewFile<-select(NewFile, patid, Diag=paste0("Diag", YearVar[i]), Any=paste0("Any", YearVar[i]))
  
  NewFile<-NewFile%>%
    group_by(Diag)%>%
    summarise(ActCount=n(), AnyCount=sum(Any), Year=as.character(paste0(YearVar[i])))  
  AllDiag<-rbind(AllDiag, NewFile)
} 

AllDiag$Diag<-as.factor(AllDiag$Diag)

levels(AllDiag$Diag)<-c("Bipolar disorder", "Other psychoses", "Schizophrenia")
AllDiag$Diag <- factor(AllDiag$Diag, levels = c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

AllDiag<-select(AllDiag, strat=Diag, Year, AnyCount, ActCount)
AllDiag$Group<-"DiagTV"

AllDiagCI<-BinomCI(AllDiag$AnyCount, AllDiag$ActCount)*100

AllDiag<-cbind(AllDiag, AllDiagCI)

AllDiag$Year<-case_when(AllDiag$Year=="0001" ~ "2000-2001",
                        AllDiag$Year=="0102" ~ "2001-2002",
                        AllDiag$Year=="0203" ~ "2002-2003",
                        AllDiag$Year=="0304" ~ "2003-2004",
                        AllDiag$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="SMI diagnosis")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/DiagBP.PNG", width = 10,   height = 10)

####Prescribed antipsychotics####
CPRD$FirstAnyPres<-pmin(CPRD$FirstAPDate, CPRD$FirstBPDate, CPRD$FirstEvidenceAPDate, CPRD$FirstEvidenceBPDate, na.rm=TRUE)
APFields<-select(CPRD, patid, FirstAnyPres)
All<-merge(x=All, y=APFields, by="patid", all.x=TRUE, all.y=FALSE)

AP<-All%>%
  mutate(AP0001 = case_when(FirstAnyPres<="2001-03-31" ~ 1,
                            TRUE ~ 0),
         AP0102 = case_when(FirstAnyPres<="2002-03-31"~ 1,
                            TRUE ~ 0),
         AP0203 = case_when(FirstAnyPres<="2003-03-31" ~ 1,
                            TRUE ~ 0),
         AP0304 = case_when(FirstAnyPres<="2004-03-31" ~ 1,
                            TRUE ~ 0),
         AP0405 = case_when(FirstAnyPres<="2005-03-31" ~ 1,
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
YearVar<-(c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708", "0809", "0910", "1011", "1112", "1213", "1314", "1415", "1516", "1617", "1718"))

AllAP<-subset(AP, Act0001==1)%>%
  select(patid, AP=AP0001, Any=Any0001)%>%
  group_by(AP)%>%
  summarise(ActCount=n(), AnyCount=sum(Any), Year="0001")

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

AllAP$Year<-case_when(AllAP$Year=="0001" ~ "2000-2001",
                      AllAP$Year=="0102" ~ "2001-2002",
                      AllAP$Year=="0203" ~ "2002-2003",
                      AllAP$Year=="0304" ~ "2003-2004",
                      AllAP$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Antipsychotics/mood stabilisers")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/APBP.PNG", width = 10,   height = 10)

####Create mega table####

AllTable<-rbind(AllCountry, AllSex, AllEthn, AllAge, AllExempt, AllQOF, AllDiag, AllAP)
write.csv(AllTable, file="DatamindCPRD/Outputs/BPsAll.csv")

AllAny<-All%>%
  subset(!is.na(FullIMD))%>%
  group_by(FullIMD)%>%
  select(starts_with("Any"))%>%
  summarise_all(list(sum))%>%
  pivot_longer(cols=starts_with("Any"), names_to = "AnyYear", values_to = "AnyCount")

AllActive<-All%>%
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

AllIMD$Year<-case_when(AllIMD$Year=="0001" ~ "2000-2001",
                       AllIMD$Year=="0102" ~ "2001-2002",
                       AllIMD$Year=="0203" ~ "2002-2003",
                       AllIMD$Year=="0304" ~ "2003-2004",
                       AllIMD$Year=="0405" ~ "2004-2005",
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
  theme(legend.position = c(0.2, 0.9),  legend.text = element_text(colour = "black", size = 16), legend.title = element_text(colour = "black", size = 18))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,by = 10))+
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs( x="Year", y = "Percentage of patients receiving blood pressure screening", fill="Index of multiple deprivation")+
  annotate("text", x = "2015-2016", y = 95, label = "Blood pressure", vjust = -0.5, size=8)+
  geom_vline(xintercept=11.5, linetype='dotted')+
  annotate("text", x = 11.5, y = 11, label = "Incentivised", vjust = -0.5, angle = 90, size=6)

ggsave("DatamindCPRD/Outputs/IMDBP.PNG", width = 10,   height = 10)