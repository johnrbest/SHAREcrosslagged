########################################################################################
# Topic:    Analysis of Longitudinal SHARE data - cross-lagged models
# Goal:     Use cognition and depression in cross lagged models
# Author:   John Best, PhD
# Last modified: 2021 September 21
########################################################################################


# 1. Load libraries and prepare data --------------------------------------

library(gtsummary)
library(ggpubr)
library(ggeffects)
library(lavaan)
library(semTools)
library(ggthemes)

source("Initial_processing.R")

# 2. Cross-lagged panel model specification ---------------------------------

#single group
clpm <- '
  #regressions
  cog_2 ~ a*cog_1 + b*age_1 + c*dep_1 + d*female
  cog_3 ~ a*cog_2 + b*age_2 + c*dep_2 + d*female
  cog_4 ~ a*cog_3 + b*age_3 + c*dep_3 + d*female
  cog_5 ~ a*cog_4 + b*age_4 + c*dep_4 + d*female
  cog_6 ~ a*cog_5 + b*age_5 + c*dep_5 + d*female
  
  dep_2 ~ e*cog_1 + f*age_1 + g*dep_1 + h*female
  dep_3 ~ e*cog_2 + f*age_2 + g*dep_2 + h*female
  dep_4 ~ e*cog_3 + f*age_3 + g*dep_3 + h*female
  dep_5 ~ a*cog_4 + b*age_4 + c*dep_4 + d*female
  dep_6 ~ a*cog_5 + b*age_5 + c*dep_5 + d*female

  
  #residual variances
  cog_2 ~~ i*cog_2
  cog_3 ~~ i*cog_3
  cog_4 ~~ i*cog_4
  cog_5 ~~ i*cog_5
  cog_6 ~~ i*cog_6
  
  dep_2~~ j*dep_2
  dep_3~~ j*dep_3
  dep_4~~ j*dep_4
  dep_5~~ j*dep_5
  dep_6~~ j*dep_6
  
  #(residual) covariances
  cog_1~~dep_1
  
  cog_2~~k*dep_2
  cog_3~~k*dep_3
  cog_4~~k*dep_4
  cog_5~~k*dep_5
  cog_6~~k*dep_6

'

#multi-group by age and sex
clpm_multi <- '
  #regressions
  cog_2 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_1 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_1
  cog_3 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_2 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_2 
  cog_4 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_3 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_3 
  cog_5 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_4 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_4 
  cog_6 ~ c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_5 + c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_5 

  dep_2 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_1 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_1 
  dep_3 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_2 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_2 
  dep_4 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_3 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_3
  dep_5 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_4 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_4
  dep_6 ~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_5 + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_5
  
  cog_1 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_1 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_1
  cog_2 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_2 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_2
  cog_3 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_3 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_3
  cog_4 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_4 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_4
  cog_5 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_5 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_5
  cog_6 ~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*age_6 + c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*chronic_mod_6
  
  dep_1 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_1 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_1
  dep_2 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_2 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_2
  dep_3 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_3 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_3
  dep_4 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_4 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_4
  dep_5 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_5 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_5
  dep_6 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_6 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_6
  
  cog_1 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_2 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_3 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_4 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_5 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_6 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  
  dep_1 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_2 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_3 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_4 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_5 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_6 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  
  #residual variances
  cog_2 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*cog_2
  cog_3 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*cog_3
  cog_4 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*cog_4
  cog_5 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*cog_5
  cog_6 ~~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*cog_6
  
  dep_2~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*dep_2
  dep_3~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*dep_3
  dep_4~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*dep_4
  dep_5~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*dep_5
  dep_6~~ c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*dep_6
  
  #(residual) covariances
  cog_1~~dep_1
  
  cog_2~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_2
  cog_3~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_3
  cog_4~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_4
  cog_5~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_5
  cog_6~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_6
'

#remove last age group
clpm_multi_under85 <- gsub(",[[:lower:]]10|,[[:lower:]]9","",clpm_multi)
clpm_multi_under85_v2 <- gsub(",[[:lower:]]10|,[[:lower:]]9|,[[:lower:]]8|,[[:lower:]]7|,[[:lower:]]6|,[[:lower:]]5|",
                              "",clpm_multi)

# 3. Random intercepts cross-lagged panel model specification---------------

#single group
riclpm <- '
  #random intercept factors
  cog =~ 1*cog_1 + 1*cog_2 + 1*cog_3 + 1*cog_4 + 1*cog_5 + 1*cog_6
  dep =~ 1*dep_1 + 1*dep_2 + 1*dep_3 + 1*dep_4 + 1*dep_5 + 1*dep_6

  cog_1~~0*cog_1
  cog_2~~0*cog_2
  cog_3~~0*cog_3
  cog_4~~0*cog_4
  cog_5~~0*cog_5
  cog_6~~0*cog_6
  
  dep_1~~0*dep_1
  dep_2~~0*dep_2
  dep_3~~0*dep_3
  dep_4~~0*dep_4
  dep_5~~0*dep_5
  dep_6~~0*dep_6
  
  #wave-specific residual factors
  cog_1r=~1*cog_1
  cog_2r=~1*cog_2
  cog_3r=~1*cog_3
  cog_4r=~1*cog_4
  cog_5r=~1*cog_5
  cog_6r=~1*cog_6
  
  dep_1r=~1*dep_1
  dep_2r=~1*dep_2
  dep_3r=~1*dep_3
  dep_4r=~1*dep_4
  dep_5r=~1*dep_5
  dep_6r=~1*dep_6
  
  cog_2r~~a*cog_2r
  cog_3r~~a*cog_3r
  cog_4r~~a*cog_4r
  cog_5r~~a*cog_5r
  cog_6r~~a*cog_6r
  
  dep_2r~~b*dep_2r
  dep_3r~~b*dep_3r
  dep_4r~~b*dep_4r
  dep_5r~~b*dep_5r
  dep_6r~~b*dep_6r
  
  #structural relations between residual factors
  
  cog_2r~ c*cog_1r + d*dep_1r
  cog_3r~ c*cog_2r + d*dep_2r
  cog_4r~ c*cog_3r + d*dep_3r
  cog_5r~ c*cog_4r + d*dep_4r
  cog_6r~ c*cog_5r + d*dep_5r

  dep_2r~ e*dep_1r + f*cog_1r
  dep_3r~ e*dep_2r + f*cog_2r
  dep_4r~ e*dep_3r + f*cog_3r
  dep_5r~ e*dep_4r + f*cog_4r
  dep_6r~ e*dep_5r + f*cog_5r
  
  cog_1r~~dep_1r
  
  cog_2r~~g*dep_2r
  cog_3r~~g*dep_3r
  cog_4r~~g*dep_4r
  cog_5r~~g*dep_5r
  cog_6r~~g*dep_6r
  
  cog + dep ~~ 0*cog_1r + 0*dep_1r

'

#multi group by age and sex
riclpm_multi <- '
  #adjust outcomes for age and other covariates
  cog_1 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_1 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_1
  cog_2 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_2 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_2
  cog_3 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_3 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_3
  cog_4 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_4 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_4
  cog_5 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_5 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_5
  cog_6 ~ c(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10)*age_6 + c(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)*chronic_mod_6
  
  dep_1 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_1 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_1
  dep_2 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_2 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_2
  dep_3 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_3 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_3
  dep_4 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_4 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_4
  dep_5 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_5 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_5
  dep_6 ~ c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)*age_6 + c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)*chronic_mod_6
  
  cog_1 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_2 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_3 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_4 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_5 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  cog_6 ~ c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10)*edu_1 + c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)*income_1
  
  dep_1 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_2 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_3 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_4 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_5 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  dep_6 ~ c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)*edu_1 + c(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)*income_1
  
  #random intercept factors
  cog =~ 1*cog_1 + 1*cog_2 + 1*cog_3 + 1*cog_4 + 1*cog_5 + 1*cog_6
  dep =~ 1*dep_1 + 1*dep_2 + 1*dep_3 + 1*dep_4 + 1*dep_5 + 1*dep_6
  
  cog_1~~0*cog_1
  cog_2~~0*cog_2
  cog_3~~0*cog_3
  cog_4~~0*cog_4
  cog_5~~0*cog_5
  cog_6~~0*cog_6
  
  dep_1~~0*dep_1
  dep_2~~0*dep_2
  dep_3~~0*dep_3
  dep_4~~0*dep_4
  dep_5~~0*dep_5
  dep_6~~0*dep_6
  
  dep ~~ c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)*cog
  
  #wave-specific residual factors
  cog_1r=~1*cog_1
  cog_2r=~1*cog_2
  cog_3r=~1*cog_3
  cog_4r=~1*cog_4
  cog_5r=~1*cog_5
  cog_6r=~1*cog_6
  
  dep_1r=~1*dep_1
  dep_2r=~1*dep_2
  dep_3r=~1*dep_3
  dep_4r=~1*dep_4
  dep_5r=~1*dep_5
  dep_6r=~1*dep_6
  
  cog_1r~~cog_1r
  dep_1r~~dep_1r
  
  cog_2r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_2r
  cog_3r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_3r
  cog_4r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_4r
  cog_5r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_5r
  cog_6r~~c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)*cog_6r
  
  dep_2r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_2r
  dep_3r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_3r
  dep_4r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_4r
  dep_5r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_5r
  dep_6r~~c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)*dep_6r
  
  #structural relations between residual factors
  
  cog_2r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_1r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_1r
  cog_3r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_2r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_2r
  cog_4r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_3r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_3r
  cog_5r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_4r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_4r
  cog_6r~ c(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)*cog_5r + c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)*dep_5r

  dep_2r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*dep_1r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*cog_1r
  dep_3r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*dep_2r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*cog_2r
  dep_4r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*dep_3r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*cog_3r
  dep_5r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*dep_4r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*cog_4r
  dep_6r~ c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)*dep_5r + c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)*cog_5r
  
  cog_1r~~dep_1r
  
  cog_2r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_2r
  cog_3r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_3r
  cog_4r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_4r
  cog_5r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_5r
  cog_6r~~c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)*dep_6r
  
  cog + dep ~~ 0*cog_1r + 0*dep_1r
'
#remove last age group
riclpm_multi_under85 <- gsub(",[[:lower:]]10|,[[:lower:]]9","",riclpm_multi)
riclpm_multi_under85_v2 <- gsub(",[[:lower:]]10|,[[:lower:]]9|,[[:lower:]]8|,[[:lower:]]7|,[[:lower:]]6|,[[:lower:]]5|",
                                "",riclpm_multi)


# 4.  Create cross-lagged functions ---------------------------------------

cross_lagged_function <- function(cogvar,depvar,cogname){
  temp_wide <- df %>% 
    as_tibble() %>% 
    rename(cog = cogvar,
           dep = depvar) %>% 
    mutate(income = thinc_m/10000) %>% 
    select(mergeid,sequence,wave,age,female,datapoints,
           #primary variables of interest
           cog,dep,
           #other covariates
           bmi,VigPA,maxgrip,chronic_mod,edu,income
    ) %>% 
    #drop_na(cog,dep) %>% 
    pivot_wider(id_cols = c(mergeid,female,datapoints,
    ),
    names_from=sequence,
    values_from=c(age,dep,cog,
                  bmi,VigPA,maxgrip,chronic_mod,edu,income)) %>% 
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
                                  `5_1`="85 to 95 female"),
           agebin_female_v2 = fct_collapse(agebin_female,
                                           "45 to 65 male" = c("45 to 55 male","55 to 65 male"),
                                           "45 to 65 female" = c("45 to 55 female","55 to 65 female"),
                                           "65+ male" = c("65 to 75 male","75 to 85 male","85 to 95 male"),
                                           "65+ female" = c("65 to 75 female","75 to 85 female","85 to 95 female"))
           )
  
  #fixed effects model
  clpm_multi_fit <- temp_wide %>% 
    filter(datapoints>2,
           #agebin<5 #exclude 85+
           ) %>% 
    arrange(agebin_female_v2) %>% 
    cfa(clpm_multi_under85_v2,
        data=.,
        missing="ML",                
        fixed.x=F,                
        group="agebin_female_v2"
    )
  cog_dep_clpm_fitmeasures[[cogvar]] <<- fitmeasures(clpm_multi_fit)
  #print(fitmeasures(clpm_multi_fit))
  cog_dep_clpm_fitsummary[[cogvar]] <<- summary(clpm_multi_fit,estimates=F)
  #print(summary(clpm_multi_fit,estimates=F))
  standardizedSolution(clpm_multi_fit) %>%
    as_tibble() %>%
    write.csv(.,file=paste0(cogvar,depvar,"_clpm_estimates.csv"))
  
  plot_clpm <- standardizedSolution(clpm_multi_fit) %>% 
    as_tibble() %>% 
    mutate(label = parameterestimates(clpm_multi_fit)[,"label"]) %>% 
    filter(str_detect(label,"b|c|g")) %>% 
    mutate(group_names = recode(group,
                                `2`="45_to_65 Male",
                                `4`="65+ Male",
                                #`6`="65_to_75 Male",
                                #`8`="75_to_85 Male",
                                #`10`="85_to_95 Male",
                                `1`="45_to_65 Female",
                                `3`="65+ Female",
                                #`5`="65_to_75 Female",
                                #`7`="75_to_85 Female",
                                #`9`="85_to_95 Female"
    )) %>% 
    mutate(path = str_replace_all(label,
                                  c("^b\\d{1,2}" = "cognition on prior depression",
                                    "^c\\d{1,2}" = "depression on prior cognition",
                                    "^g\\d{1,2}" = "residual covariation"
                                  ))
    ) %>% 
    group_by(group,path) %>% 
    filter(row_number()==1) %>% 
    separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
    mutate(Age = str_replace_all(Age,"_"," "),
           path = factor(path,
                         levels = c("depression on prior cognition",
                                    "cognition on prior depression",
                                    "residual covariation"
                         ),
                         labels = c("Path a: prior cognition -> depression",
                                    "Path b: prior depression -> cognition",
                                    "Path c: depression <-> cognition"
                         ))
    ) %>% 
    ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std,colour=Age,linetype=Sex)) +
    geom_point(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper),position=position_dodge(0.3),width=.2) +
    theme_classic(base_size=10) +
    facet_wrap(~path,scales = "free") +
    ylab("Estimate") +
    ylim(-.2,0) +
    xlab("") +
    geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_flip() +
    scale_colour_colorblind() +
    ggtitle(paste0("Depression and ",cogname))
  
  print(plot_clpm)
  ggsave(plot_clpm,file=paste0("clpm_",cogvar,"_dep_plot.png"),height=7,width=10)
  
  #random intercepts model
  riclpm_multi_fit <- temp_wide %>% 
    filter(datapoints>3,
           #agebin<5
           ) %>% 
    arrange(agebin_female_v2) %>% 
    cfa(riclpm_multi_under85_v2,
        data=.,
        missing="ML",                
        fixed.x=F,
        group="agebin_female_v2")
  
  cog_dep_riclpm_fitmeasures[[cogvar]] <<- fitmeasures(riclpm_multi_fit)
  #print(fitmeasures(riclpm_multi_fit))
  cog_dep_riclpm_fitsummary[[cogvar]] <<- summary(riclpm_multi_fit,estimates=F)
  #print(summary(riclpm_multi_fit,estimates=F))
  standardizedSolution(riclpm_multi_fit) %>%
    as_tibble() %>%
    write.csv(.,file=paste0(cogvar,depvar,"_riclpm_estimates.csv"))
  
  plot_riclpm <- standardizedSolution(riclpm_multi_fit) %>% 
    as_tibble() %>%
    mutate(label = parameterestimates(riclpm_multi_fit)[,"label"]) %>%
    filter(str_detect(label,"d|f|g|h")) %>% 
    mutate(group_names = recode(group,
                                `2`="45_to_65 Male",
                                `4`="65+ Male",
                                #`6`="65_to_75 Male",
                                #`8`="75_to_85 Male",
                                #`10`="85_to_95 Male",
                                `1`="45_to_65 Female",
                                `3`="65+ Female",
                                #`5`="65_to_75 Female",
                                #`7`="75_to_85 Female",
                                #`9`="85_to_95 Female"
    )) %>% 
    mutate(path = str_replace_all(label,
                                  c("^d\\d{1,2}" = "cognition on prior depression (within-person)",
                                    "^f\\d{1,2}" = "depression on prior cognition (within-person)",
                                    "^g\\d{1,2}" = "residual covariation (within-person)",
                                    "^h\\d{1,2}" = "covariation (between-person)"))
    ) %>% 
    group_by(group,path) %>% 
    filter(row_number()==1) %>% 
    separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
    mutate(Age = str_replace_all(Age,"_"," "),
           path = factor(path,
                         levels = c("depression on prior cognition (within-person)",
                                    "cognition on prior depression (within-person)",
                                    "residual covariation (within-person)",
                                    "covariation (between-person)"),
                         labels = c("Path a: prior cognition -> depression",
                                    "Path b: prior depression -> cognition",
                                    "Path c: depression <-> cognition (Within)",
                                    "Path d: depression <-> cognition (Between)"))
    ) %>% 
    ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std,colour=Age,linetype=Sex)) +
    geom_point(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper),position=position_dodge(0.3),width=.2) +
    theme_classic(base_size=10) +
    facet_wrap(~path,scales = "free") +
    ylab("Estimate") +
    ylim(-.4,.02) +
    xlab("") +
    geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_flip() +
    scale_colour_colorblind() +
    ggtitle(paste0("Depression and ",cogname))
  
  print(plot_riclpm)
  ggsave(plot_riclpm,file=paste0("riclpm_",cogvar,"_dep_plot.png"),height=7,width=10)
}


# 5. Descriptives ------------------------------------------------------------

#first assessment
df %>% 
  filter(age>=45,
         age<95,
         datapoints>3,
         sequence==1) %>%
  select(age,female,
         #datapoints,
         thinc_m,
         edu,
         chronic_mod,
         eurod,
         recall_1,recall_2,animal_fluency) %>% 
  tbl_summary(
    label = list(female ~ "Gender, female",
                 #datapoints ~ "Number of assessments completed",
                 chronic_mod ~ "Number of chronic diseases",
                 edu ~ "Educational attainment",
                 recall_1 ~ "Immediate recall: number of words",
                 recall_2 ~ "Delayed recall: number of words"
    ),
    missing_text = "(Missing)",
    type = list(chronic_mod ~ "continuous2",
                all_continuous() ~ "continuous2"),
    statistic = all_continuous() ~ c("{median} [{p25}, {p75}]",
                                     "{mean} ({sd})")
  ) %>%
  modify_header(label = "**Variable**") %>% 
  as_flex_table()

#time-varying
df %>% 
  filter(age>=45,
         age<95,
         datapoints>3,
         sequence<7) %>% 
  select(sequence,
         eurod,
         recall_1,recall_2,
         animal_fluency) %>% 
  tbl_summary(
    by = sequence,
    label = list(
                 recall_1 ~ "Immediate recall: number of words",
                 recall_2 ~ "Delayed recall: number of words"
                 ),
    missing_text = "(Missing)",
    type = list(all_continuous() ~ "continuous2"),
    statistic = all_continuous() ~ c("{median} [{p25}, {p75}]",
                                     "{mean} ({sd})")
  ) %>%
  modify_header(label = "**Variable**") %>% 
  modify_spanning_header(all_stat_cols() ~ "**Assessment wave**") %>% 
  as_flex_table()

# 6. Run cross-lagged functions ----------------------------------------------

cog_dep_riclpm_fitmeasures <- cog_dep_riclpm_fitsummary <- list()
cog_dep_clpm_fitmeasures <- cog_dep_clpm_fitsummary <- list()

cross_lagged_function("recall_1","logeurod","Immediate recall")
cross_lagged_function("recall_2","logeurod","Delayed recall")
cross_lagged_function("animal_fluency","logeurod","Animal fluency")

write.csv(cog_dep_clpm_fitmeasures,file="cog_dep_clpm_fitmeasures.csv")
write.csv(cog_dep_riclpm_fitmeasures,file="cog_dep_riclpm_fitmeasures.csv")

# 7.  Use lmer to quantify stability --------------------------------------

library(broom.mixed)
library(lmerTest)

outcomes <- c("logeurod","recall_1","recall_2","animal_fluency")
for (outcome in outcomes){
  mdl0 <- df %>% 
    filter(age>=45,
           age<95,
           datapoints>3,
           sequence<7) %>% 
    mutate(
      agebin=cut(age,breaks=seq(from=45,to=95,by=10),right=FALSE,labels=FALSE),
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
                                  `5_1`="85 to 95 female"),
           agebin_female_v2 = fct_collapse(agebin_female,
                                           "45 to 65 male" = c("45 to 55 male","55 to 65 male"),
                                           "45 to 65 female" = c("45 to 55 female","55 to 65 female"),
                                           "65+ male" = c("65 to 75 male","75 to 85 male","85 to 95 male"),
                                           "65+ female" = c("65 to 75 female","75 to 85 female","85 to 95 female"))
    ) %>% 
    group_by(agebin_female_v2) %>% 
    do(mdl_icc = performance::icc(lmer(as.formula(paste0(outcome,"~ 1 + (1|mergeid)")),data=.))
       )
  print(mdl0$mdl_icc)
}

for (outcome in outcomes){
  mdl0 <- df %>% 
    filter(age>=45,
           age<95,
           datapoints>3,
           sequence<7) %>% 
    lmer(as.formula(paste0(outcome,"~ 1 + (1|mergeid)")),data=.)
  print(performance::icc(mdl0))
}


# 8. Create combined plots across cognitive measures ----------------------

#CLPM
recall_1logeurod_clpm_estimates <- read_csv("recall_1logeurod_clpm_estimates.csv") %>% 
  filter(str_detect(label,"b|c|g")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^b\\d{1,2}" = "cognition on prior depression",
                                  "^c\\d{1,2}" = "depression on prior cognition",
                                  "^g\\d{1,2}" = "residual covariation"
                                )),
         Cognition = "Immediate recall"
  ) 

recall_2logeurod_clpm_estimates <- read_csv("recall_2logeurod_clpm_estimates.csv") %>% 
  filter(str_detect(label,"b|c|g")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^b\\d{1,2}" = "cognition on prior depression",
                                  "^c\\d{1,2}" = "depression on prior cognition",
                                  "^g\\d{1,2}" = "residual covariation"
                                )),
         Cognition = "Delayed recall"
  ) 

animal_fluencylogeurod_clpm_estimates <- read_csv("animal_fluencylogeurod_clpm_estimates.csv") %>% 
  filter(str_detect(label,"b|c|g")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^b\\d{1,2}" = "cognition on prior depression",
                                  "^c\\d{1,2}" = "depression on prior cognition",
                                  "^g\\d{1,2}" = "residual covariation"
                                )),
         Cognition = "Animal fluency"
  )   

combined_plot_clpm <- bind_rows(recall_1logeurod_clpm_estimates,recall_2logeurod_clpm_estimates) %>% 
  bind_rows(.,animal_fluencylogeurod_clpm_estimates) %>% 
  group_by(group,path,Cognition) %>% 
  filter(row_number()==1) %>% 
  separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
  mutate(Age = str_replace_all(Age,"_"," "),
         path = factor(path,
                       levels = c("depression on prior cognition",
                                  "cognition on prior depression",
                                  "residual covariation"
                       ),
                       labels = c("Path a: prior cognition -> depression",
                                  "Path b: prior depression -> cognition",
                                  "Path c: depression <-> cognition"
                       ))
  ) %>% 
  arrange(Age,Sex,path) %>% 
  ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std)) +
  geom_point(aes(shape=Cognition,colour=Age),position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper,colour=Age,linetype=Sex),position=position_dodge(0.3),width=.2) +
  theme_classic(base_size=10) +
  facet_wrap(~path,scales = "free") +
  ylab("Standardized Estimate") +
  ylim(-.2,0) +
  xlab("") +
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  scale_colour_colorblind()
ggsave(combined_plot_clpm,file="combined_plot_clpm.png",height=6.5,width=10)

#RI-CLPM
recall_1logeurod_riclpm_estimates <- read_csv("recall_1logeurod_riclpm_estimates.csv") %>% 
  filter(str_detect(label,"d|f|g|h")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^d\\d{1,2}" = "cognition on prior depression (within-person)",
                                  "^f\\d{1,2}" = "depression on prior cognition (within-person)",
                                  "^g\\d{1,2}" = "residual covariation (within-person)",
                                  "^h\\d{1,2}" = "covariation (between-person)")),
         Cognition = "Immediate recall"
  )

recall_2logeurod_riclpm_estimates <- read_csv("recall_2logeurod_riclpm_estimates.csv") %>% 
  filter(str_detect(label,"d|f|g|h")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^d\\d{1,2}" = "cognition on prior depression (within-person)",
                                  "^f\\d{1,2}" = "depression on prior cognition (within-person)",
                                  "^g\\d{1,2}" = "residual covariation (within-person)",
                                  "^h\\d{1,2}" = "covariation (between-person)")),
         Cognition = "Delayed recall"
  )

animal_fluencylogeurod_riclpm_estimates <- read_csv("animal_fluencylogeurod_riclpm_estimates.csv") %>% 
  filter(str_detect(label,"d|f|g|h")) %>% 
  select(-c(1)) %>% 
  mutate(group_names = recode(group,
                              `2`="45_to_65 Male",
                              `4`="65+ Male",
                              #`6`="65_to_75 Male",
                              #`8`="75_to_85 Male",
                              #`10`="85_to_95 Male",
                              `1`="45_to_65 Female",
                              `3`="65+ Female",
                              #`5`="65_to_75 Female",
                              #`7`="75_to_85 Female",
                              #`9`="85_to_95 Female"
  )) %>% 
  mutate(path = str_replace_all(label,
                                c("^d\\d{1,2}" = "cognition on prior depression (within-person)",
                                  "^f\\d{1,2}" = "depression on prior cognition (within-person)",
                                  "^g\\d{1,2}" = "residual covariation (within-person)",
                                  "^h\\d{1,2}" = "covariation (between-person)")),
         Cognition = "Animal fluency"
  )

combined_plot_riclpm <- bind_rows(recall_1logeurod_riclpm_estimates,recall_2logeurod_riclpm_estimates) %>% 
  bind_rows(.,animal_fluencylogeurod_riclpm_estimates) %>% 
  group_by(group,path,Cognition) %>%
  filter(row_number()==1) %>% 
  separate(group_names, into = c("Age","Sex"),sep = " ") %>% 
  mutate(Age = str_replace_all(Age,"_"," "),
         path = factor(path,
                       levels = c("depression on prior cognition (within-person)",
                                  "cognition on prior depression (within-person)",
                                  "residual covariation (within-person)",
                                  "covariation (between-person)"),
                       labels = c("Path a: prior cognition -> depression",
                                  "Path b: prior depression -> cognition",
                                  "Path c: depression <-> cognition (Within)",
                                  "Path d: depression <-> cognition (Between)"))
  ) %>% 
  arrange(Age,Sex,path) %>% 
  ggplot(aes(x = order(Age,group,decreasing=TRUE),y=est.std,colour=Age,linetype=Sex,shape=Cognition,)) +
  geom_point(position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin = ci.lower,ymax = ci.upper),position=position_dodge(0.3),width=.2) +
  theme_classic(base_size=10) +
  facet_wrap(~path,scales = "free",nrow = 1) +
  ylab("Standardized Estimate") +
  ylim(-.4,.02) +
  xlab("") +
  geom_hline(yintercept=0,linetype="dashed",alpha=0.5) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  scale_colour_colorblind()
ggsave(combined_plot_riclpm,file="combined_plot_riclpm.png",height=6.5,width=12)
