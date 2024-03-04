## PREP SUBJECTS AND PHENOTYPES FOR UKB ## 

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(data.table)
library(readxl);library(naniar)

#Import connectome subject id - this is after QC
ukb_conn_subj <- read.csv("UKBiobank_8185_roi_vol.csv", header = T)[1]
colnames(ukb_conn_subj) <- "connectome_id"

#import linkage data 
linkage_id <- read.csv("bridge_file_10279_4844_GALE_MCINTOSH.csv", header = F)
colnames(linkage_id) <- c("connectome_id","f.eid") 

#import mdd cidi
mdd_cidi <- readRDS("MHQ.1907.ukb24262.Derived.rds") %>% 
  select(f.eid, Depressed.Ever) %>% drop_na()
colnames(mdd_cidi)[2] <- "mdd_cidi" 

#to exclude neurological disorders and psychiatric disorders following earlier publication
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5514104/
exclusions <- readRDS("ukb24262-mdd.mdd_exclusions.rds")
colnames(exclusions)[2:ncol(exclusions)] <- paste0("excl_",colnames(exclusions)[2:ncol(exclusions)])
keep <- exclusions %>% 
  filter(excl_parkinsons_icd == F & excl_bipolar_icd == F & excl_mpd_icd == F & 
           excl_scz_icd == F & excl_autism_icd == F & excl_id_icd == F) %>% 
  select(f.eid)

rm_connid <- c("1497301", "4928353", "2394167", "2934482") #withdrawn/disconnected

#determine final sample
to_keep <- merge(ukb_conn_subj, linkage_id, by = "connectome_id") %>% #N=8185
  merge(., mdd_cidi, by = "f.eid") %>% #N=5124
  .[.$f.eid %in% keep$f.eid,] %>%  #N=5107
  .[.$connectome_id %in% rm_connid == F,] %>% #N=5104
  select(connectome_id, f.eid, mdd_cidi) 
table(to_keep$mdd_cidi) #controls:3599, cases:1505

#save
saveRDS(to_keep,"ukb_mdd_pheno_compiled.rds")


#Prep covariates ==== 

#2021
baseline <- readRDS("BaselineCharacteristics.rds")
recruitment <- readRDS("Recruitment.rds")
touchscreen <- readRDS("Touchscreen.rds")
physical <- readRDS("PhysicalMeasures.rds")

#sex
sex <- baseline %>% .[,c(1,2)]
colnames(sex)[2] <- "sex"
sex$sex <- ifelse(sex$sex == "Male", 1,0)
sex$sex <- as.factor(sex$sex)

#age and site
#age: age at instance 2 (age during visit to img assessment centre)
#site: 4 centres for imaging assessment but this sample came from one site 
#Cheadle (11025), Reading (11026), Newcastle (11027), Bristol (11028)
mri_age_site <- recruitment %>% select(f.eid, f.54.2.0,f.21003.2.0)
colnames(mri_age_site)[2:3]<- c("site.id","mri_age")
mri_age_site$site.id <- as.factor(mri_age_site$site.id)

#ethnicity 
ethnic <- touchscreen %>% select(f.eid, f.21000.0.0)
colnames(ethnic)[2] <- "ethnicity"

#smoking status
smoking <- touchscreen %>% select(f.eid,f.20116.0.0, f.20116.2.0) #take instance 2 to correspond to imaging visit 
colnames(smoking)[2:3] <- c("smoke_baseline","smoke_img")
smoking <- smoking %>% naniar::replace_with_na(replace = list(smoke_baseline = c("Prefer not to answer"),
                                                              smoke_img = c("Prefer not to answer")))
smoking$smoke_baseline <- gsub("Never",0, smoking$smoke_baseline)
smoking$smoke_baseline <- gsub("Previous",0, smoking$smoke_baseline) 
smoking$smoke_baseline <- gsub("Current",1, smoking$smoke_baseline)
smoking$smoke_img <- gsub("Never",0, smoking$smoke_img)
smoking$smoke_img <- gsub("Previous",0, smoking$smoke_img)
smoking$smoke_img <- gsub("Current",1, smoking$smoke_img)
# apply(smoking[2:3], 2, function(x) table(x))

smoking$smoke_baseline <- as.factor(smoking$smoke_baseline)
smoking$smoke_img <- as.factor(smoking$smoke_img)

#alcohol status
alcohol <- touchscreen %>% select(f.eid, f.20117.0.0,f.20117.2.0)
colnames(alcohol)[2:3] <- c("alcohol_baseline","alcohol_img")
alcohol <- alcohol %>% naniar::replace_with_na(replace = list(alcohol_baseline = c("Prefer not to answer"),
                                                              alcohol_img = c("Prefer not to answer")))
alcohol$alcohol_baseline <- gsub("Never",0, alcohol$alcohol_baseline)
alcohol$alcohol_baseline <- gsub("Previous",0, alcohol$alcohol_baseline) #follow STRADL and put previous as 0?
alcohol$alcohol_baseline <- gsub("Current",1, alcohol$alcohol_baseline)
alcohol$alcohol_img <- gsub("Never",0, alcohol$alcohol_img)
alcohol$alcohol_img <- gsub("Previous",0, alcohol$alcohol_img)
alcohol$alcohol_img <- gsub("Current",1, alcohol$alcohol_img)
# apply(alcohol[2:3], 2, function(x) table(x))

alcohol$alcohol_baseline <- as.factor(alcohol$alcohol_baseline)
alcohol$alcohol_img <- as.factor(alcohol$alcohol_img)

#bmi
bmi <- physical %>% select(f.eid, f.21001.0.0, f.21001.2.0)
colnames(bmi)[2:3] <- c("bmi_baseline","bmi_img")

#education
edu <- touchscreen %>% select(f.eid, contains("f.6138.0")) 

edu[-1] <- apply(edu[-1], 2, function(col){
  stringi::stri_replace_all_regex(col,
                                  pattern=c('^College or University degree$', 
                                            '^A levels\\/AS levels or equivalent$',
                                            '^O levels\\/GCSEs or equivalent$', 
                                            '^CSEs or equivalent$',
                                            '^NVQ or HND or HNC or equivalent$', 
                                            '^Other professional qualifications eg: nursing, teaching$',
                                            '^None of the above$',
                                            '^Prefer not to answer$'),
                                  replacement=c(6,5,4,3,2,1,NA,NA),
                                  vectorize=FALSE)
}) #works better than nested gsub

edu$highest_edu <- apply(edu[-1], 1, max, na.rm = T) #select the highest value for each subject
edu$highest_edu_bin <- ifelse(edu$highest_edu == 7, 1, 0) #dichotomize to w/wo university degree
edu <- edu %>% select(f.eid, highest_edu, highest_edu_bin)

#income
income <- touchscreen %>% select(f.eid,f.738.0.0,f.738.2.0) 
income[-1] <- apply(income[-1], 2, function(col){
  stringi::stri_replace_all_regex(col,
                                  pattern=c('^Greater than 100,000$', 
                                            '^52,000 to 100,000$',
                                            '^31,000 to 51,999$', 
                                            '^18,000 to 30,999$',
                                            '^Less than 18,000$', 
                                            '^Do not know$',
                                            '^Prefer not to answer$'),
                                  replacement=c(5,4,3,2,1,NA,NA),
                                  vectorize=FALSE)
}) 
colnames(income) <- c("f.eid", "income_baseline", "income_img")

#antidepressant status (file prepped by Mark)
antidep <- read_csv("UKB_MHQMDD.csv") %>% select(f.eid, DrugDepression)
colnames(antidep)[2] <- "antidep"

#ICV 
#extracted from field 26521-2.0 "Volume of EstimatedTotalIntraCranial (whole brain)" 
#on Eddie cuz file size too big to do locally
icv <- readRDS("ukb_icv.rds") 

#Merge all covariates

#full sample
ukb_covs <- full_join(sex, mri_age_site) %>% 
  full_join(., ethnic) %>% 
  full_join(., smoking) %>% 
  full_join(., alcohol) %>% 
  full_join(., bmi) %>% 
  full_join(., edu) %>% 
  full_join(., income) %>% 
  full_join(., antidep) %>% 
  full_join(., icv) 

#conn_only sample 
ukb_covs_conn <- merge(ukb_covs, to_keep[, c("connectome_id", "f.eid")], by = "f.eid")

#save RDS 
saveRDS(ukb_covs_conn,"ukb_covs_conn_only.rds")

