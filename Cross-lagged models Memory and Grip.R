########################################################################################
# Topic:    Analysis of Longitudinal SHARE data
# Author:   John Best, PhD
# Last modified: 2021 February 05
########################################################################################


# 1. Load libraries and prepare data --------------------------------------

# Need to first execute easySHARE_rel7_0_0.rda file

library(tidyverse)
library(rms)
library(tableone)
library(ggpubr)
library(ggeffects)
library(lavaan)
library(semTools)
library(ggthemes)

load("easySHARE_rel7_1_0.rda")

data <- easySHARE_rel7_1_0
data[data<0] <- NA
data <- data %>% 
  mutate(retired=ifelse(ep005_>1,0,1))

mean_age_data <- data %>% 
  group_by(mergeid) %>% 
  summarise(mean_age=mean(age,na.rm=TRUE),
            min_age=min(age),
            max_age=max(age),
            datapoints=sum(!is.na(age)),
            sex=mean(female))

essentials <- c("mergeid","coupleid","wave","int_year","int_month","age","female",
                "bmi","retired","thinc_m","isced1997_r","mar_stat","chronic_mod")
cog <- c("recall_1","recall_2") #primary cognitive variable of interest: delayed recall
cov <- c("mobilityind","eurod","maxgrip","br015_","adla","iadlza") #cov refers to the primary correlate varaible to examine alongside cognition


#Transformations and norming
df <- data[,c(essentials,cog,cov)] %>%
  arrange(wave) %>% 
  mutate(cubedmobilityind=mobilityind^(1/3),
         logeurod=log(eurod+1),
         cubedadla=adla^(1/3),
         cubediadlza=iadlza^(1/3),
         VigPA=5-br015_,
         Sex=ifelse(female==0,"Male","Female"),
         agebin=cut(age,breaks=seq(from=45,to=95,by=2),right=FALSE,labels=FALSE),
         number=1,
         edu = isced1997_r
  ) %>%
  group_by(agebin,Sex) %>%
  mutate(recall_2_norm=as.numeric(scale(recall_2)),
         cubedmobilityind_norm=as.numeric(scale(cubedmobilityind)),
         logeurod_norm=as.numeric(scale(logeurod)),
         VigPA_norm=as.numeric(scale(VigPA)),
         maxgrip_norm=as.numeric(scale(maxgrip)),
         bmi_norm=as.numeric(scale(bmi))
  ) %>% 
  ungroup() %>% 
  group_by(mergeid) %>% 
  mutate(timepoint=cumsum(number),
         datapoints=sum(!is.na(age)),
         firstvisit=min(wave,na.rm=TRUE),
         lastvisit=max(wave,na.rm=TRUE))

df$longitudinal <- ifelse(df$datapoints>1,"Yes","No") #Distinguish those with 2+ timepoints from those without
df$edu[df$edu==95|df$edu==97] <- NA #Remove "other" and "still in school" from education variable

#create variable labels
label(df[["age"]]) <- "Age"
label(df[["edu"]]) <- "Edu"
label(df[["int_year"]]) <- "Year"
label(df[["wave"]]) <- "Wave"
label(df[["chronic_mod"]]) <- "Chronic diseases"
label(df[["maxgrip_norm"]]) <- "Grip"
label(df[["VigPA_norm"]]) <- "PA"
label(df[["cubedmobilityind_norm"]]) <- "Mobility"
label(df[["logeurod_norm"]]) <- "Depression"
label(df[["bmi_norm"]]) <- "BMI"

dd <- datadist(df);options(datadist='dd')

#descriptives for first visit
vars <- c("wave","int_year","datapoints","age","retired","thinc_m","isced1997_r","mar_stat","chronic_mod")
catvars <- c("wave","int_year","retired","isced1997_r","mar_stat")
df[,vars] <- apply(df[,vars],2,as.numeric)
print(CreateTableOne(vars=vars,strata="Sex",factorVars=catvars,data=subset(df,wave==firstvisit&datapoints>1),test=FALSE,addOverall=TRUE),
      quote=TRUE,nonnormal=c("thinc_m","chronic_mod","datapoints"),smd=TRUE)

#descriptives for last visit (longitudinal only)
print(CreateTableOne(vars=vars,strata="Sex",factorVars=catvars,data=subset(df,wave==lastvisit&datapoints>1),test=FALSE,addOverall=TRUE),
      quote=TRUE,nonnormal=c("thinc_m","chronic_mod"),smd=TRUE)

#compare those with longitudinal data to those without
vars <- c("age","Sex","retired","thinc_m","isced1997_r","mar_stat","chronic_mod","recall_2","bmi","eurod","VigPA","maxgrip")
catvars <- c("Sex","retired","isced1997_r","mar_stat")
print(CreateTableOne(vars=vars,strata=c("longitudinal"),factorVars=catvars,data=subset(df,wave==firstvisit),
                     test=FALSE),quote=TRUE,nonnormal=c("thinc_m","chronic_mod"),smd=TRUE)


# 2. Transform to long format ---------------------------------------------

#Use grip strength and delayed memory recall

df_wide <- df %>% 
  as_tibble() %>% 
  select(mergeid,timepoint,wave,age,female,recall_2,maxgrip,datapoints) %>% 
  drop_na() %>% #drop timepoints with missing data on above
  rename(recall2=recall_2) %>% 
  pivot_wider(id_cols = c(mergeid,female,datapoints),
              names_from=timepoint,
              values_from=c(age,maxgrip,recall2)) %>% 
  filter(age_1>=45,
         age_1<95) %>% 
  mutate(
    agebin=cut(age_1,breaks=seq(from=45,to=95,by=10),right=FALSE,labels=FALSE),
    agebin_female = paste(agebin,female,sep="_")
  ) %>% 
  mutate(agebin_female = recode(agebin_female,
                                `1_0`="45 to 55 male",
                                `2_0`="55 to 65 male",
                                `3_0`="65 to 75 male",
                                `4_0`="75 to 85 male",
                                `5_0`="85 to 95 male",
                                `1_1`="45 to 55 female",
                                `2_1`="55 to 65 female",
                                `3_1`="65 to 75 female",
                                `4_1`="75 to 85 female",
                                `5_1`="85 to 95 female"
                                ))

# 3. Cross-lagged panel model ---------------------------------------------

#single group
clpm <- '
  #regressions
  recall2_2 ~ a*recall2_1 + b*age_1 + c*maxgrip_1 + d*female
  recall2_3 ~ a*recall2_2 + b*age_2 + c*maxgrip_2 + d*female
  recall2_4 ~ a*recall2_3 + b*age_3 + c*maxgrip_3 + d*female
  recall2_5 ~ a*recall2_4 + b*age_4 + c*maxgrip_4 + d*female
  
  maxgrip_2 ~ e*recall2_1 + f*age_1 + g*maxgrip_1 + h*female
  maxgrip_3 ~ e*recall2_2 + f*age_2 + g*maxgrip_2 + h*female
  maxgrip_4 ~ e*recall2_3 + f*age_3 + g*maxgrip_3 + h*female
  maxgrip_5 ~ e*recall2_4 + f*age_4 + g*maxgrip_4 + h*female

  
  #residual variances
  recall2_2 ~~ i*recall2_2
  recall2_3 ~~ i*recall2_3
  recall2_4 ~~ i*recall2_4
  recall2_5 ~~ i*recall2_5
  
  maxgrip_2~~ j*maxgrip_2
  maxgrip_3~~ j*maxgrip_3
  maxgrip_4~~ j*maxgrip_4
  maxgrip_5~~ j*maxgrip_5
  
  #(residual) covariances
  recall2_1~~maxgrip_1
  
  recall2_2~~k*maxgrip_2
  recall2_3~~k*maxgrip_3
  recall2_4~~k*maxgrip_4
  recall2_5~~k*maxgrip_5

'

clpm_fit <- df_wide %>% 
  filter(datapoints>=2) %>% 
  cfa(clpm,
      data=.,
      missing="ML",
      fixed.x=FALSE)

#fitmeasures(clpm_fit)
#summary(clpm_fit)
#standardizedSolution(clpm_fit)

#multi-group by age and sex
clpm_multi <- '
  #regressions
  recall2_2 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_1 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_1
  recall2_3 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_2 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_2 
  recall2_4 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_3 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_3 
  recall2_5 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_4 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_4 

  maxgrip_2 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_1 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_1 
  maxgrip_3 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_2 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_2 
  maxgrip_4 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_3 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_3
  maxgrip_5 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_4 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_4
  
  recall2_1 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_1
  recall2_2 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_2
  recall2_3 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_3
  recall2_4 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_4
  recall2_5 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_5
  
  maxgrip_1 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_1
  maxgrip_2 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_2
  maxgrip_3 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_3
  maxgrip_4 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_4
  maxgrip_5 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_5

  #residual variances
  recall2_2 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*recall2_2
  recall2_3 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*recall2_3
  recall2_4 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*recall2_4
  recall2_5 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*recall2_5
  
  maxgrip_2~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*maxgrip_2
  maxgrip_3~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*maxgrip_3
  maxgrip_4~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*maxgrip_4
  maxgrip_5~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*maxgrip_5
  
  #(residual) covariances
  recall2_1~~maxgrip_1
  
  recall2_2~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_2
  recall2_3~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_3
  recall2_4~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_4
  recall2_5~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_5
'
clpm_multi_under85 <- gsub(",[[:lower:]]10|,[[:lower:]]9","",clpm_multi)

clpm_multi_fit <- df_wide %>% 
  filter(datapoints>=2,
         agebin<5) %>% 
  arrange(agebin_female) %>% 
  cfa(clpm_multi_under85,
      data=.,
      missing="ML",                
      fixed.x=FALSE,                
      group="agebin_female")

fitmeasures(clpm_multi_fit)
summary(clpm_multi_fit,estimates=F)

clpm_memory_grip_plot <- standardizedSolution(clpm_multi_fit) %>% 
  as_tibble() %>% 
  mutate(label = parameterestimates(clpm_multi_fit)[,"label"]) %>% 
  filter(str_detect(label,"b|c|g")) %>% 
  mutate(group_names = recode(group,
                                `2`="45_to_55 Male",
                                `4`="55_to_65 Male",
                                `6`="65_to_75 Male",
                                `8`="75_to_85 Male",
                                `10`="85_to_95 Male",
                                `1`="45_to_55 Female",
                                `3`="55_to_65 Female",
                                `5`="65_to_75 Female",
                                `7`="75_to_85 Female",
                                `9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^b\\d{1,2}" = "recall on prior grip",
                                  "^c\\d{1,2}" = "grip on prior recall",
                                  "^g\\d{1,2}" = "residual covariation"
                                  ))
  ) %>% 
  group_by(group,path) %>% 
  filter(row_number()==1) %>% 
  separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
  mutate(Age = str_replace_all(Age,"_"," "),
         path = factor(path,
                       levels = c("grip on prior recall",
                                  "recall on prior grip",
                                  "residual covariation"
                                  ),
                       labels = c("Prior memory -> Grip",
                                  "Prior grip -> Memory",
                                  "Grip <-> Memory"
                                  ))
  ) %>% 
  ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std,colour=Age,linetype=Sex)) +
  geom_point(position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper),position=position_dodge(0.3),width=.2) +
  theme_classic(base_size=10) +
  facet_wrap(~path,scales = "free") +
  ylab("Estimate") +
  ylim(-.2,.3) +
  xlab("") +
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  scale_colour_colorblind()
ggsave(clpm_memory_grip_plot,file="clpm_memory_grip_plot.png",height=7,width=10)

# 4. Random intercepts cross-lagged panel model ---------------------------

#single group
riclpm <- '
  #random intercept factors
  recall2 =~ 1*recall2_1 + 1*recall2_2 + 1*recall2_3 + 1*recall2_4 + 1*recall2_5
  maxgrip =~ 1*maxgrip_1 + 1*maxgrip_2 + 1*maxgrip_3 + 1*maxgrip_4 + 1*maxgrip_5

  recall2_1~~0*recall2_1
  recall2_2~~0*recall2_2
  recall2_3~~0*recall2_3
  recall2_4~~0*recall2_4
  recall2_5~~0*recall2_5
  
  maxgrip_1~~0*maxgrip_1
  maxgrip_2~~0*maxgrip_2
  maxgrip_3~~0*maxgrip_3
  maxgrip_4~~0*maxgrip_4
  maxgrip_5~~0*maxgrip_5
  
  #wave-specific residual factors
  recall2_1r=~1*recall2_1
  recall2_2r=~1*recall2_2
  recall2_3r=~1*recall2_3
  recall2_4r=~1*recall2_4
  recall2_5r=~1*recall2_5
  
  maxgrip_1r=~1*maxgrip_1
  maxgrip_2r=~1*maxgrip_2
  maxgrip_3r=~1*maxgrip_3
  maxgrip_4r=~1*maxgrip_4
  maxgrip_5r=~1*maxgrip_5
  
  recall2_2r~~a*recall2_2r
  recall2_3r~~a*recall2_3r
  recall2_4r~~a*recall2_4r
  recall2_5r~~a*recall2_5r
  
  maxgrip_2r~~b*maxgrip_2r
  maxgrip_3r~~b*maxgrip_3r
  maxgrip_4r~~b*maxgrip_4r
  maxgrip_5r~~b*maxgrip_5r
  
  #structural relations between residual factors
  
  recall2_2r~ c*recall2_1r + d*maxgrip_1r
  recall2_3r~ c*recall2_2r + d*maxgrip_2r
  recall2_4r~ c*recall2_3r + d*maxgrip_3r
  recall2_5r~ c*recall2_4r + d*maxgrip_4r

  maxgrip_2r~ e*maxgrip_1r + f*recall2_1r
  maxgrip_3r~ e*maxgrip_2r + f*recall2_2r
  maxgrip_4r~ e*maxgrip_3r + f*recall2_3r
  maxgrip_5r~ e*maxgrip_4r + f*recall2_4r
  
  recall2_1r~~maxgrip_1r
  
  recall2_2r~~g*maxgrip_2r
  recall2_3r~~g*maxgrip_3r
  recall2_4r~~g*maxgrip_4r
  recall2_5r~~g*maxgrip_5r
  
  recall2 + maxgrip ~~ 0*recall2_1r + 0*maxgrip_1r

'

riclpm_fit <- df_wide %>% 
  filter(datapoints>=2) %>% 
  cfa(riclpm,
      data=.,
      missing="ML",
      fixed.x=FALSE
      )

#fitmeasures(riclpm_fit)
#summary(riclpm_fit)
#standardizedSolution(riclpm_fit)

#multi group by age and sex
riclpm_multi <- '
  #adjust outcomes for age
  recall2_1 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_1
  recall2_2 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_2
  recall2_3 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_3
  recall2_4 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_4
  recall2_5 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_5
  
  maxgrip_1 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_1
  maxgrip_2 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_2
  maxgrip_3 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_3
  maxgrip_4 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_4
  maxgrip_5 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_5
  
  #random intercept factors
  recall2 =~ 1*recall2_1 + 1*recall2_2 + 1*recall2_3 + 1*recall2_4 + 1*recall2_5
  maxgrip =~ 1*maxgrip_1 + 1*maxgrip_2 + 1*maxgrip_3 + 1*maxgrip_4 + 1*maxgrip_5
  
  recall2_1~~0*recall2_1
  recall2_2~~0*recall2_2
  recall2_3~~0*recall2_3
  recall2_4~~0*recall2_4
  recall2_5~~0*recall2_5
  
  maxgrip_1~~0*maxgrip_1
  maxgrip_2~~0*maxgrip_2
  maxgrip_3~~0*maxgrip_3
  maxgrip_4~~0*maxgrip_4
  maxgrip_5~~0*maxgrip_5
  
  maxgrip ~~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*recall2
  
  #wave-specific residual factors
  recall2_1r=~1*recall2_1
  recall2_2r=~1*recall2_2
  recall2_3r=~1*recall2_3
  recall2_4r=~1*recall2_4
  recall2_5r=~1*recall2_5
  
  maxgrip_1r=~1*maxgrip_1
  maxgrip_2r=~1*maxgrip_2
  maxgrip_3r=~1*maxgrip_3
  maxgrip_4r=~1*maxgrip_4
  maxgrip_5r=~1*maxgrip_5
  
  recall2_1r~~recall2_1r
  maxgrip_1r~~maxgrip_1r
  
  recall2_2r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_2r
  recall2_3r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_3r
  recall2_4r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_4r
  recall2_5r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*recall2_5r
  
  maxgrip_2r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_2r
  maxgrip_3r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_3r
  maxgrip_4r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_4r
  maxgrip_5r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*maxgrip_5r
  
  #structural relations between residual factors
  
  recall2_2r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_1r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_1r
  recall2_3r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_2r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_2r
  recall2_4r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_3r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_3r
  recall2_5r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*recall2_4r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*maxgrip_4r

  maxgrip_2r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*maxgrip_1r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*recall2_1r
  maxgrip_3r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*maxgrip_2r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*recall2_2r
  maxgrip_4r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*maxgrip_3r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*recall2_3r
  maxgrip_5r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*maxgrip_4r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*recall2_4r
  
  recall2_1r~~maxgrip_1r
  
  recall2_2r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_2r
  recall2_3r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_3r
  recall2_4r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_4r
  recall2_5r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*maxgrip_5r
  
  recall2 + maxgrip ~~ 0*recall2_1r + 0*maxgrip_1r
'
riclpm_multi_under85 <- gsub(",[[:lower:]]10|,[[:lower:]]9","",riclpm_multi)

riclpm_multi_fit <- df_wide %>% 
  filter(datapoints>=2,
         agebin<5) %>% 
  arrange(agebin_female) %>% 
  cfa(riclpm_multi_under85,
      data=.,
      missing="ML",                
      fixed.x=FALSE,
      group="agebin_female")

fitmeasures(riclpm_multi_fit)
summary(riclpm_multi_fit,estimates=F)


riclpm_memory_grip_plot <- standardizedSolution(riclpm_multi_fit) %>% 
  as_tibble() %>%
  mutate(label = parameterestimates(riclpm_multi_fit)[,"label"]) %>%
  filter(str_detect(label,"d|f|g|h")) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_55 Male",
                              `4`="55_to_65 Male",
                              `6`="65_to_75 Male",
                              `8`="75_to_85 Male",
                              `10`="85_to_95 Male",
                              `1`="45_to_55 Female",
                              `3`="55_to_65 Female",
                              `5`="65_to_75 Female",
                              `7`="75_to_85 Female",
                              `9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                            c("^d\\d{1,2}" = "recall on prior grip (within-person)",
                              "^f\\d{1,2}" = "grip on prior recall (within-person)",
                              "^g\\d{1,2}" = "residual covariation (within-person)",
                              "^h\\d{1,2}" = "covariation (between-person)"))
         ) %>% 
  group_by(group,path) %>% 
  filter(row_number()==1) %>% 
  separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
  mutate(Age = str_replace_all(Age,"_"," "),
         path = factor(path,
                       levels = c("grip on prior recall (within-person)",
                                  "recall on prior grip (within-person)",
                                  "residual covariation (within-person)",
                                  "covariation (between-person)"),
                       labels = c("Prior memory -> Grip",
                                  "Prior grip -> Memory",
                                  "Grip <-> Memory (Within)",
                                  "Grip <-> Memory (Between)"))
         ) %>% 
  ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std,colour=Age,linetype=Sex)) +
  geom_point(position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper),position=position_dodge(0.3),width=.2) +
  theme_classic(base_size=10) +
  facet_wrap(~path,scales = "free") +
  ylab("Estimate") +
  ylim(-.3,.7) +
  xlab("") +
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  scale_colour_colorblind()
ggsave(riclpm_memory_grip_plot,file="riclpm_memory_grip_plot.png",height=7,width=10)

