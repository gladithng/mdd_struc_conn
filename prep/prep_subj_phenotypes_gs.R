## PREP SUBJECTS AND PHENOTYPES FOR GS ## 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(data.table)
library(readxl);library(naniar)

#Import connectome subject id - this is after QC
stradl_conn_subj <- read.csv("subjects.csv", header = F) #N=938 
colnames(stradl_conn_subj) <- "ID"

#import linkage data
linkage_id <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")

#import mdd cidi
#GS recontact
load("STRADL.RData") #compiled previously; should see x in env
mdd_cidi <- x %>%  
  select(id,CIDI_MDD) %>% #7667 controls, 1506 cases
  filter(CIDI_MDD == 0 | CIDI_MDD == 1) %>% #Control = 0, MDD = 1
  merge(linkage_id[2:3], ., by.x = "genotype_ID", by.y = "id")  %>% 
  select(stradl_ID, CIDI_MDD) #736 controls, 168 cases
colnames(mdd_cidi) <- c("ID","cidi_GS")

#exclusions - follow ukb
basic_info <- read_excel("STRADL_collated.xlsx", sheet = "basicinfo")
neuro <- basic_info[grep("parkinson's", basic_info$IllnessDiagnosis, ignore.case = T),] 
psych <- basic_info[grep("bipolar|multiple personality|schizophrenia|autism|intellectual disability", basic_info$IllnessDiagnosis, ignore.case = T),] 

rm <- basic_info[1] %>% 
  mutate(., to_exclude = ifelse(.$ID %in% neuro$ID | .$ID %in% psych$ID, 1,0))
table(rm$to_exclude)

#determine final sample 
to_keep <- merge(stradl_conn_subj, linkage_id, by.x = "ID", by.y = "stradl_ID") %>% #N=938
  merge(., mdd_cidi, by = "ID") %>% #N=729
  merge(., rm, by = "ID") %>% 
  filter(to_exclude == 0) %>% #N=725
  select(ID, cidi_GS)
table(to_keep$cidi_GS) #controls: 593, cases: 132

#save
saveRDS(to_keep,"stradl_mdd_pheno_compiled.rds")

#Prep covariates ====
#R_:recontact assessment;P_:pre-appointment assessment;C_:clinical assessment

#demographic variables at imaging assessment
demo <- read_excel("STRADL_collated.xlsx", sheet = "demographics") %>% 
  select(ID,AgeFaceToFace,Sex)
colnames(demo) <- c("ID","age","sex")

#site
site <- read_excel("ST id linkage final.xlsx") %>% select(st,centre)
site$centre <- as.factor(site$centre)
colnames(site) <- c("ID", "site")

#lifestyle factors - smoking and alcohol
alcohol <- read_excel("STRADL_collated.xlsx", sheet = "alcohol") %>% 
  select(ID, CurrentDrinker)
colnames(alcohol)[2] <- "alcohol_current"
alcohol$alcohol_current <- as.factor(alcohol$alcohol_current)

smoke <- read_excel("STRADL_collated.xlsx", sheet = "tobacco") %>% 
  select(ID, EverSmoker, CurrentSmoker) %>% 
  mutate_at(vars(EverSmoker:CurrentSmoker), as.factor)
colnames(smoke) <- c("ID", "smoke_ever", "smoke_current")

#icv (file by Aleks)
icv <- read_csv("STRADL_Measures_Standardised_ICV.csv") 
colnames(icv)[2] <- "icv_std"

#antidepressants - follow those in UKB
med <- read_excel("STRADL_collated.xlsx", sheet = "basicinfo") %>% 
  select(ID, CurrentsMeds) %>% 
  mutate(mlt = strsplit(as.character(CurrentsMeds), ",")) %>% #split the commas into separate rows
  unnest(mlt)

ukb_antidep <- read.delim("ls.antidepressants_ukb_code.txt")
id <- med$ID[tolower(med$mlt) %in% tolower(ukb_antidep$antidep)] %>% unique()
antidep <- data.frame(ID = unique(med$ID))
antidep$antidep <- ifelse(antidep$ID %in% id, 1, 0)
antidep$antidep <- as.factor(antidep$antidep)

#bmi 
physical <- read_excel("STRADL_collated.xlsx", sheet = "PhysicalMeasures") %>% 
  select(ID, `Height(cm)`, `Weight(kg)`)
colnames(physical) <- c("ID", "height", "weight")
physical$bmi <- physical$weight / (physical$height/100)^2
bmi <- physical %>% select(ID, bmi)

#education 
edu <- read_excel("STRADL_collated.xlsx",sheet = "basicinfo") %>% 
  select(ID,Education) %>% mutate_at(vars(Education),as.factor)
edu$edu_bin <- ifelse(edu$Education == 3, 1, 0) #group post secondary into 1 group
colnames(edu) <- c("ID","edu", "edu_bin")

#income (from genscot timepoint)
load("masterDB.RData") #should see totaldata in env

gs_income <- totaldata %>% 
  select(id, income) %>% #(1 - less than 10,000, 2 - between 10,000 and 30,000, 3 - between 30,000 and 50,000, 4 - between 50,000 and 70,000, 5 - more than 70,000, 6 - prefer not to answer, 99 - not known)
  mutate_at(vars(income), as.factor) %>% 
  merge(., linkage_id[2:3], by.x="id", by.y="genotype_ID") %>% 
  select(stradl_ID, income)
colnames(gs_income) <- c("ID","GS_income")

gs_income$GS_income <- as.factor(gsub("6","NA",gs_income$GS_income))
gs_income$GS_income <- as.factor(gsub("7","NA",gs_income$GS_income))
gs_income$GS_income <- as.factor(gsub("99","NA",gs_income$GS_income))
table(gs_income$GS_income)

#ethnicity
ethnic <- totaldata %>% select(id, race) %>% 
  merge(linkage.id[,c(2:3)],., by.x="genotype_ID", by.y="id") %>% 
  select(stradl_ID,race)
colnames(ethnic) <- c("ID","ethnicity")

#Merge all covariates

#full sample
stradl_covs <- full_join(demo, site) %>% 
  full_join(., ethnic) %>% 
  full_join(., smoke) %>% 
  full_join(., alcohol) %>% 
  full_join(., bmi) %>% 
  full_join(., edu) %>% 
  full_join(., gs_income) %>% 
  full_join(., antidep) %>% 
  full_join(., icv)

#conn only sample
stradl_covs_conn <- stradl_covs[stradl_covs$ID %in% to_keep$ID,]

saveRDS(stradl_covs_conn,"stradl_covs_conn_only.rds")

