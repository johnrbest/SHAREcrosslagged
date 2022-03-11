########################################################################################
# Topic:    Analysis of Longitudinal SHARE data - cross-lagged models
# Goal:     Conduct initial data processing
# Author:   John Best, PhD
# Last modified: 2021 September 21
########################################################################################


# Load libraries and prepare data --------------------------------------

library(tidyverse)
library(haven)

# Need to first execute easySHARE_rel7_0_0.rda file

load("easySHARE_rel7_1_0.rda")

#Import relevant Wave 8 data
health_w8 <- read_sav("sharew8_rel1-0-0_gv_health.sav") %>% 
  select(mergeid,bmi,cf008tot,cf016tot,
         numeracy,numeracy2,
         orienti,
         chronicw8c,eurod,iadl,
         maxgrip) %>% 
  rename(recall_1 = cf008tot,
         recall_2 = cf016tot,
         numeracy_1 = numeracy,
         numeracy_2 = numeracy2,
         chronic_mod = chronicw8c
  )
isced_w8 <- read_sav("sharew8_rel1-0-0_gv_isced.sav") %>% 
  select(mergeid,isced1997_r)
retired_data_w8 <- read_sav("sharew8_rel1-0-0_ep.sav") %>% 
  select(mergeid,ep005_)
demos_w8 <- read_sav("sharew8_rel1-0-0_dn.sav") %>% 
  select(mergeid,dn014_) %>% 
  rename(mar_stat = dn014_)
cover_screen_w8 <- read_sav("sharew8_rel1-0-0_cv_r.sav") %>% 
  select(mergeid,coupleid8,
         int_year,int_month,country,age_int,gender) %>% 
  mutate(wave = 8,
         female = ifelse(gender==2,1,
                         ifelse(gender==1,0,NA))) %>% 
  rename(age = age_int,
         coupleid = coupleid8)
behavioral_risks_w8 <- read_sav("sharew8_rel1-0-0_br.sav") %>% 
  select(mergeid,br015_)
cognition_w8 <- read_sav("sharew8_rel1-0-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluecy = cf010_)

wave8_data <- cover_screen_w8 %>% 
  full_join(.,health_w8,by="mergeid") %>% 
  full_join(.,behavioral_risks_w8,by="mergeid") %>% 
  full_join(.,isced_w8,by="mergeid") %>% 
  full_join(.,retired_data_w8,by="mergeid") %>%
  full_join(.,demos_w8,by="mergeid") %>% 
  full_join(.,cognition_w8,by="mergeid")

#important animal fluency data from previous timepoints

cognition_w1 <- read_sav("sharew1_rel7-1-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 1)
cognition_w2 <- read_sav("sharew2_rel7-1-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 2)
cognition_w4 <- read_sav("sharew4_rel7-1-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 4)
cognition_w5 <- read_sav("sharew5_rel7-1-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 5)
cognition_w6 <- read_sav("sharew6_rel7-1-0_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 6)
cognition_w7 <- read_sav("sharew7_rel7-1-1_cf.sav") %>% 
  select(mergeid,cf010_) %>% 
  rename(animal_fluency = cf010_) %>% 
  mutate(wave = 7)

cognition_data <- cognition_w1 %>% 
  bind_rows(cognition_w2) %>% 
  bind_rows(cognition_w4) %>%
  bind_rows(cognition_w5) %>%
  bind_rows(cognition_w6) %>%
  bind_rows(cognition_w7) 

# combine data

data <- easySHARE_rel7_1_0 %>% 
  full_join(cognition_data,by=c("mergeid","wave")) %>% 
  bind_rows(wave8_data)

data[data<0] <- NA
data <- data %>% 
  mutate(retired=ifelse(ep005_>1,0,1))

# Data transformations ----------------------------------------------------

essentials <- c("mergeid","coupleid","wave","int_year","int_month","age","female",
                "bmi","retired","thinc_m","isced1997_r","mar_stat","chronic_mod")
cog <- c("recall_1","recall_2","animal_fluency") #primary cognitive variable of interest: delayed recall
cov <- c("mobilityind","eurod","maxgrip","br015_","adla","iadlza") #cov refers to the primary correlate varaible to examine alongside cognition

df <- data[,c(essentials,cog,cov)] %>%
  group_by(mergeid) %>% 
  mutate(year = as.numeric(int_year),
         prior_time = lag(year,order_by=wave),
         datapoints=sum(!is.na(year)),
         firstvisit=min(wave,na.rm=TRUE),
         lastvisit=max(wave,na.rm=TRUE),
         sequence = row_number()
  ) %>% 
  ungroup() %>% 
  mutate(cubedmobilityind=mobilityind^(1/3),
         logeurod=log(eurod+1),
         cubedadla=adla^(1/3),
         cubediadlza=iadlza^(1/3),
         VigPA=5-br015_,
         Sex=ifelse(female==0,"Male","Female"),
         agebin=cut(age,breaks=seq(from=45,to=95,by=2),right=FALSE,labels=FALSE),
         number=1,
         edu = isced1997_r
  ) 

df$edu[df$edu==95|df$edu==97] <- NA #Remove "other" and "still in school" from education variable

