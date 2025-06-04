rm(list = ls(all.names = TRUE))
####Libraries####
library(haven)
library(psych)
library(tidyverse)
library(tableone)
library(stringr)

####Load SMI cohort####

load("DatamindCPRD/Data/CleanCPRD.Rdata")

####Load BP lists####
BPQOF<-read.table("CodeLists/Blood Pressure/BPQOF.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".")
BPAurum<-read.table("CodeLists/Blood Pressure/BPAurum.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".", colClasses=c(MedCodeId="character"))

length(which(BPQOF$Read_Code %in% BPAurum$CleansedReadCode))
BPQOF$Code_Description<-tolower(BPQOF$Code_Description)
BPAurum$Term<-tolower(BPAurum$Term)

BPQOFCodes<-merge(x=BPAurum, y=BPQOF, by.x="Term", by.y="Code_Description", all.x=FALSE, all.y=TRUE)
BPQOFCodes<-subset(BPQOFCodes, CleansedReadCode!=""| is.na(CleansedReadCode))
BPQOFCodes$MedCodeId[BPQOFCodes$Read_Code=="246A."]<-"619931000006119"
BPQOFCodes$MedCodeId[BPQOFCodes$Read_Code=="2469."]<-"114311000006111"

BPAurum$QOF<-0
BPAurum$QOF[BPAurum$MedCodeId %in% BPQOFCodes$MedCodeId]<-1

BPGold<-read.table("CodeLists/Blood Pressure/BPGold.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".", colClasses=c(medcode="character"))

length(which(BPQOF$Read_Code %in% BPGold$readcode))

BPGold$readterm<-tolower(BPGold$readterm)

BPQOFCodesG<-merge(x=BPGold, y=BPQOF, by.x="readterm", by.y="Code_Description", all.x=FALSE, all.y=TRUE)
BPGold$QOF<-0
BPGold$QOF[BPGold$medcode %in% BPQOFCodesG$medcode]<-1

CheckA<-subset(BPAurum, !(Term %in% BPGold$readterm) & CleansedReadCode!="")
CheckG<-subset(BPGold, !(readterm %in% BPAurum$Term))

BPGold<-subset(BPGold, readcode!="9OD3.00")
BPAurum<-subset(BPAurum, OriginalReadCode!="8I3Y")
BPAurum<-subset(BPAurum, OriginalReadCode!="8OAH")

GoldBPExtra<-read.table("Codelists/Blood pressure/BPGoldExtra.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".", colClasses=c(medcode="character"))
GoldBPExtra$QOF<-0
BPGold<-rbind(BPGold, GoldBPExtra)

write.table(BPGold, file = "Codelists/Blood pressure/BPGoldFinal.txt", col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")
write.table(BPAurum, file = "Codelists/Blood pressure/BPAurumFinal.txt", col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")

####Add in new codes and determine whether screening or not####
#aurum
AurumCodes<-read.table("CodeLists/Blood pressure/BPAurumFinal.txt", header=TRUE, fill=TRUE, colClasses=c(MedCodeId="character"))
AurumCodes$type<-"Screening"

AurumMed<-read_dta("Dictionary Browsers/CPRD_CodeBrowser_201904_Aurum/CPRDAurumMedical.dta")
AurumMed$term<-tolower(AurumMed$term)
AurumMed$term<-gsub('"', '', AurumMed$term)

AurumMed<-subset(AurumMed, !(originalreadcode %in% AurumCodes$OriginalReadCode))

#Note: this was the first attempt at generating code lists in R - I now have better work flow including documentation of exclusion terms!

AurumMed<-AurumMed%>%
  mutate(BP = case_when(grepl("bp|blood pressure|hypoten|hypoten|systolic|diastolic", term) ~ "BP",
                        TRUE ~ "drop"))

NewCodes<-subset(AurumMed, BP == "BP")
NewCodes<-subset(NewCodes, medcodeid=="3135013" | medcodeid=="53452019" | medcodeid=="93494011" | medcodeid=="128722012" | medcodeid=="99042012" | 
                   medcodeid=="75071013" | medcodeid=="253532015" | medcodeid=="264473010" | medcodeid=="264472017" | medcodeid=="80224019" | 
                   medcodeid=="264471012" | medcodeid=="99047018" | medcodeid=="131046010" | medcodeid=="284062018" | medcodeid=="300835018" | 
                   medcodeid=="251674014" | medcodeid=="299678018" | medcodeid=="299684015" | medcodeid=="285673019" | medcodeid=="285657013" | medcodeid=="285662014" | 
                   medcodeid=="350841017" | medcodeid=="405053013" | medcodeid=="2548266018" | medcodeid=="443764015" | medcodeid=="1780252014" | 
                   medcodeid=="1780318013" | medcodeid=="285265015" |  medcodeid=="1780253016" | 
                   medcodeid=="1780319017" |medcodeid=="411919011" | medcodeid=="350601000000116" | 
                   medcodeid=="884961000006115" | medcodeid=="909441000006118" | medcodeid=="1131851000000118" | medcodeid=="790031000006119" | medcodeid=="1576591000006118" | medcodeid=="1760931000000113" | 
                   medcodeid=="1806081000006115" | medcodeid=="1806071000006118" | medcodeid=="1877391000006112" | medcodeid=="1858771000006114" | medcodeid=="2193971000000110" | 
                   medcodeid=="1806141000006113" | medcodeid=="1855781000006115" | medcodeid=="1916581000006119" | medcodeid=="1916631000006116" | medcodeid=="1978701000006112" | medcodeid=="2193031000000112" | medcodeid=="2193021000000110" | 
                   medcodeid=="2115801000000110" | medcodeid=="1858801000006111" | medcodeid=="2403261000000112" | medcodeid=="1916641000006114" | medcodeid=="1916621000006119" | medcodeid=="1992471000006116" | medcodeid=="2117971000000119" | medcodeid=="2403221000000116" |
                   medcodeid=="2403141000000118" |medcodeid=="940001000006110" |medcodeid=="299677011" |medcodeid=="300871017" |medcodeid=="300870016" |
                   medcodeid=="299680012" |medcodeid=="299682016" |medcodeid=="300988013" |medcodeid=="395751018" |medcodeid=="213111000006112" |
                   medcodeid=="151161000006115" |medcodeid=="158241000006117" |medcodeid=="351361000000117" |medcodeid=="153941000006118" |medcodeid=="498791000006113" |
                   medcodeid=="351341000000118" |medcodeid=="790091000006115" |medcodeid=="790121000006116" |medcodeid=="790141000006111" |medcodeid=="790071000006116" |medcodeid=="790131000006118" |
                   medcodeid=="1846991000006111" |medcodeid=="264487010" |medcodeid=="1846941000006119" |medcodeid=="264485019" |medcodeid=="2474335018" |
                   medcodeid=="1846961000006115" |medcodeid=="164121000006118" |medcodeid=="285658015" |medcodeid=="299681011" |medcodeid=="1138321000000113" |
                   medcodeid=="523921000006115" |medcodeid=="1706551000000114" |medcodeid=="1823901000006112" |medcodeid=="300837014" |medcodeid=="461309011" |medcodeid=="884121000006111" |
                   medcodeid=="285661019" |medcodeid=="351421000000117" | medcodeid=="262960017" | medcodeid=="19421000006119" | medcodeid=="19411000006110" | medcodeid=="19421000006119" | 
                   medcodeid=="19431000006116" | medcodeid=="19441000006114" | medcodeid=="19451000006111" | medcodeid=="2423811000000118"
                 | medcodeid=="305764016" | medcodeid=="212961000006118" | medcodeid=="305818015"|medcodeid=="504901000006118"| medcodeid=="504911000006115"| medcodeid=="1409014" 
                 | medcodeid=="300869017" | medcodeid=="451424017"  | medcodeid=="109700019"  | medcodeid=="60444016"  | medcodeid=="299654011"  
                 | medcodeid=="264467014" | medcodeid=="789941000006117"  | medcodeid=="789951000006115"  | medcodeid=="789961000006118"  | medcodeid=="789971000006113" 
                 | medcodeid=="789991000006114"  | medcodeid=="790001000006110"  | medcodeid=="741701000006114"  | medcodeid=="790011000006113"  | medcodeid=="395753015" 
                 | medcodeid=="1855791000006117"  | medcodeid=="856971000006112" | medcodeid=="84112010"  | medcodeid=="64168014"  | medcodeid=="299687010"
                 | medcodeid=="299687010"  | medcodeid=="143003017"  | medcodeid=="299675015"  | medcodeid=="299673010"  | medcodeid=="107545013"  
                 | medcodeid=="299655012"  | medcodeid=="741661000006118" | medcodeid=="741661000006118"  | medcodeid=="741691000006114"  | medcodeid=="1908721000006111"|  
                   medcodeid=="884131000006114"|medcodeid=="64282015"| medcodeid=="299665018"| medcodeid=="264486018"| medcodeid=="2160047013"| 
                   medcodeid=="110659019"| medcodeid=="90135019"| medcodeid=="299650019"| medcodeid=="728671000006119"| medcodeid=="728671000006119"| 
                   medcodeid=="108730018"| medcodeid=="213081000006118"| medcodeid=="213091000006115"| medcodeid=="213101000006114"| medcodeid=="887811000006117"|
                   medcodeid=="1908711000006115"| medcodeid=="530161000000111"| medcodeid=="535651000000113"| medcodeid=="530221000000112"| medcodeid=="535591000000118"|
                   medcodeid=="2159168015"| medcodeid=="741681000006111" | medcodeid=="299686018"| 
                   medcodeid=="15881000000111" | medcodeid=="15861000000119" | medcodeid=="21601000000111" | medcodeid=="16051000000119" | 
                   medcodeid=="15841000000115" | medcodeid=="15901000000114" | medcodeid=="26091000000116")
#whats left
Check<-subset(AurumMed, !(medcodeid %in% NewCodes$medcodeid) & BP=="BP")

NewCodes<-select(NewCodes, -BP)
NewCodes$QOF<-0
NewCodes<-NewCodes%>%
  mutate(type = case_when(medcodeid=="285265015" |medcodeid=="1760931000000113" |medcodeid=="1916581000006119" |medcodeid=="411919011" |medcodeid=="405053013" | medcodeid=="1780253016" | medcodeid=="1780319017" | medcodeid=="790031000006119" | medcodeid=="1855781000006115" | medcodeid=="2115801000000110"
                          | medcodeid=="264467014" | medcodeid=="790141000006111" |medcodeid=="790091000006115" | medcodeid=="2117971000000119" | medcodeid=="940001000006110" | medcodeid=="153941000006118" | medcodeid=="790071000006116" | medcodeid=="164121000006118" | medcodeid=="1706551000000114"| medcodeid=="461309011"  ~ "Screening",
                          TRUE ~ "value"))

NewCodes<-select(NewCodes, MedCodeId = medcodeid, Term = term, OriginalReadCode = originalreadcode, CleansedReadCode = cleansedreadcode, SnomedCTConceptId = snomedctconceptid, SnomedCTDescriptionId = snomedctdescriptionid, Release = release, QOF, type)
NaomiAurum<-rbind(AurumCodes, NewCodes)
#Remove the couple that are target
NaomiAurum<-subset(NaomiAurum, MedCodeId!="460132010" & MedCodeId!="460131015" & MedCodeId!="856881000006115" & MedCodeId!="856891000006117" )


write.table(NaomiAurum, file = "CodeLists/Blood pressure/NaomiBPAurumFinalScreenAndValue.txt",  col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")

#Gold
GoldCodes<-read.table("CodeLists/Blood pressure/BPGoldFinal.txt", header=TRUE, colClasses=c(medcode="character"))
GoldCodes$type<-"Screening"

GoldMed<-read.csv("pathfinder/Naomi/CodeLists/medicalGold.csv", header=FALSE, fill=TRUE)
GoldMed$V7<-tolower(GoldMed$V7)
GoldMed$V1<-as.character(GoldMed$V1)

GoldMed<-subset(GoldMed, !(V1 %in% GoldCodes$medcode))
GoldMed<-GoldMed%>%
  mutate(BP = case_when(grepl("bp|blood pressure|hypoten|hyperten|systolic|diastolic", V7) ~ "BP",
                        TRUE ~ "drop"))

NewCodesGold<-subset(GoldMed, BP == "BP")
NewCodesGold<-subset(NewCodesGold,  V1=="31341" | V1=="10424" | V1=="97533" | V1=="52431" | V1=="102458" | V1=="45149" | V1=="1894" | V1=="109942" | 
                       V1=="109942" | V1=="99234" | V1=="34231" | V1=="16478" | V1=="34244" | V1=="66589" | V1=="87862" | V1=="59609" | V1=="83473"| 
                       V1=="799" | V1=="10818" | V1=="16565" | V1=="2666" | V1=="85944" | V1=="103874" | V1== "99748" | V1=="96743" | V1=="102406" | V1=="19070" | 
                       V1=="12680"  | V1=="32976" | V1=="36305" | V1=="24127" | V1=="13186" | V1=="4444" | V1=="3712" | V1=="105480" | V1=="19342" | 
                       V1=="3269" | V1=="803" | V1=="31341" | V1=="34744" | V1=="109611" | V1=="110631" | V1=="18482" |V1=="22333" | 
                       V1=="30776" | V1=="3756" | V1== "32609" | V1=="32338" | V1=="98230" | V1=="15377" | V1=="18590" | V1=="3425" | V1=="72030" | 
                       V1=="62432" | V1=="73586" | V1=="27511" | V1=="93055" | V1=="43664" | V1=="44549" | V1=="107704" | V1=="5513" | V1=="4552" | V1=="57288" | 
                       V1=="51635" | V1=="7329" | V1=="42229" | V1=="105989" | V1=="105316" | V1=="105371" | V1=="105274" | V1=="22356"
                     | V1=="16059" | V1=="31755" | V1=="73293" | V1=="27634" | V1=="4344" | V1=="105487"
                     | V1=="4372" | V1=="109771" | V1=="14856" | V1=="66567" | V1=="101649"| V1=="8296"| V1=="62718" | V1=="18765" | V1=="11056")
#whats left
Check<-subset(GoldMed, !(V1 %in% NewCodesGold$V1) & BP=="BP")
CheckAurum<-subset(NewCodesGold, !(V7 %in% NewCodes$Term))
CheckGold<-subset(NewCodes, !(Term %in% NewCodesGold$V7))
CheckGold<-subset(CheckGold, CleansedReadCode!="")
NewGold<-subset(GoldMed, V2 %in% CheckGold$CleansedReadCode)
NewCodesGold<-rbind(NewCodesGold, NewGold)

NewCodesGold<-select(NewCodesGold, -BP)
NewCodesGold$QOF<-0
NewCodesGold<-NewCodesGold%>%
  mutate(type = case_when(V1=="10424" |V1=="8296" |V1=="103874" |V1=="102406"|V1=="19070"|V1=="36305"
                          |V1=="24127"|V1=="13186"|V1=="4444"|V1=="803"|V1=="109611"|V1=="110631"
                          |V1=="18482"|V1=="21826"|V1=="12948"|V1=="4552"|V1=="27634"|V1=="27634"~ "Screening",
                          TRUE ~ "value"))
GoldScreen<-subset(NewCodesGold, type=="Screening")
AurumScreen<-subset(NewCodes, type=="Screening")
CheckAurum<-subset(GoldScreen, !(V7 %in% AurumScreen$Term))
CheckGold<-subset(AurumScreen, !(Term %in% GoldScreen$V7))
CheckGold<-subset(CheckGold, CleansedReadCode!="")

NewCodesGold<-select(NewCodesGold, medcode = V1, clinicalevents = V3, referralevents=V4, testevents=V5, immunisationevents=V6, readcode=V2, readterm=V7, databasebuild = V8, QOF, type)
NaomiGold<-rbind(GoldCodes, NewCodesGold)
#Remove the couple that are a target
NaomiGold<-subset(NaomiGold, medcode!="51357" & medcode!="33330")

write.table(NaomiGold, file = "CodeLists/Blood pressure/NaomiBPGoldFinalScreenAndValue.txt",  col.names=TRUE, row.names=FALSE, sep="\t", dec = ".")
