rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(tableone)

####Load SMI cohort####

load("DatamindCPRD/Data/CleanFinalCohort.Rdata")

CPRD<-select(CPRD, -DeathTime)

#Sort region
CPRD$country<-"England"
CPRD$country[CPRD$region=="Northern Ireland"]<-"Northern Ireland"
CPRD$country[CPRD$region=="Wales"]<-"Wales"
CPRD$country[CPRD$region=="Scotland"]<-"Scotland"

CPRD$EnglandRegion<-CPRD$region
CPRD$EnglandRegion[CPRD$region=="Northern Ireland"|CPRD$region=="Wales"|CPRD$region=="Scotland"]<-NA
table(CPRD$region)
table(CPRD$country)
table(CPRD$EnglandRegion)

#Sort active period
CPRD<-CPRD%>%
  mutate(Act0001 = case_when(end>="2001-03-31" & as.Date("2001-03-31") - start>=90 ~ 1,
                             TRUE ~ 0),
    Act0102 = case_when(end>="2002-03-31" & as.Date("2002-03-31") - start>=90 ~ 1,
                        TRUE ~ 0),
    Act0203 = case_when(end>="2003-03-31" & as.Date("2003-03-31") - start>=90 ~ 1,
                        TRUE ~ 0),
    Act0304 = case_when(end>="2004-03-31" & as.Date("2004-03-31") - start>=90 ~ 1,
                        TRUE ~ 0),
    Act0405 = case_when(end>="2005-03-31" & as.Date("2005-03-31") - start>=90 ~ 1,
                             TRUE ~ 0),
         Act0506 = case_when(end>="2006-03-31" & as.Date("2006-03-31") - start>=90 ~ 1,
                             TRUE ~ 0),
         Act0607 = case_when(end>="2007-03-31" & as.Date("2007-03-31") - start>=90 ~ 1,
                             TRUE ~ 0),
         Act0708 = case_when(end>="2008-03-31" & as.Date("2008-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act0809 = case_when(end>="2009-03-31" & as.Date("2009-03-31") - start>=90  ~ 1,
                             TRUE ~ 0), 
         Act0910 = case_when(end>="2010-03-31" & as.Date("2010-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1011 = case_when(end>="2011-03-31" & as.Date("2011-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1112 = case_when(end>="2012-03-31" & as.Date("2012-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1213 = case_when(end>="2013-03-31" & as.Date("2013-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1314 = case_when(end>="2014-03-31" & as.Date("2014-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1415 = case_when(end>="2015-03-31" & as.Date("2015-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1516 = case_when(end>="2016-03-31" & as.Date("2016-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1617 = case_when(end>="2017-03-31" & as.Date("2017-03-31") - start>=90 ~ 1,
                             TRUE ~ 0), 
         Act1718 = case_when(end>="2018-03-31" & as.Date("2018-03-31") - start>=90 ~ 1,
                             TRUE ~ 0))

save(CPRD, file = "DatamindCPRD/Data/FinalCohortActive.Rdata")
