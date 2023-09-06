#Analysis that accompanies manuscript
#Title: Genetic variation in parasite avoidance, yet no evidence for constitutive fitness costs
#Authors: Caroline R. Amoroso, Leila L. Shepard, & Amanda K. Gibson

#Last updated: 6 September 2023

#relies upon the following data files in the folder "Data"
#"Lawn-leaving assay.csv"
#"Population growth assay.csv"
#"Zhang et al 2021 fecundity data.csv"

#load packages
library(tidyverse)
library(emmeans)
library(lme4)
library(car)
library(lmtest)
library(lmerTest)
library(pbkrtest)

#Read in lawn-leaving data
avoid.dat = read.csv("./Data/Lawn-leaving assay.csv", header=T)

#convert to date / time
avoid.dat$Date_time_picked = lubridate::as_datetime(avoid.dat$Date_time_picked, format="%m/%d/%y %H:%M",tz="EST") #
avoid.dat$Date_time_counted = lubridate::as_datetime(avoid.dat$Date_time_counted, format="%m/%d/%y %H:%M",tz="EST") #format="%m/%d/%y %H:%M",
avoid.dat$Time_interval = difftime(avoid.dat$Date_time_counted, avoid.dat$Date_time_picked, units="hours")

#create categorical variable of early (18h)/late (24h) count
avoid.dat$Count_time = NA
avoid.dat$Count_time[avoid.dat$Time_interval<22]="early"
avoid.dat$Count_time[avoid.dat$Time_interval>=22]="late"
avoid.dat$Count_time = factor(avoid.dat$Count_time)

#Summary stats on the amount of time between picking and counting
mean(avoid.dat$Time_interval[avoid.dat$Count_time=="early"])
sd(avoid.dat$Time_interval[avoid.dat$Count_time=="early"])

mean(avoid.dat$Time_interval[avoid.dat$Count_time=="late"])
sd(avoid.dat$Time_interval[avoid.dat$Count_time=="late"])

#set baseline factor for Ce strain (N2), bacterial strain (OP50)
avoid.dat$C.elegans_strain = factor(avoid.dat$C.elegans_strain, levels=c("N2", "CX11314", "JU258", "MY23", "CB4856",
                                                                         "EG4725", "ED3017", "DL238", "JT11398", 
                                                                         "LKC34", "MY16", "JU775"))
avoid.dat$Bacterial_strain = factor(avoid.dat$Bacterial_strain, levels=c("OP50","Db10"))

#how many counts did we not account for all 20 worms?
length(which(avoid.dat$Worms_ON+avoid.dat$Worms_OFF+avoid.dat$Excl_from_sub20 <20)) / length(avoid.dat$Worms_ON)
length(which(avoid.dat$Worms_ON+avoid.dat$Worms_OFF+avoid.dat$Excl_from_sub20 <20 & avoid.dat$Bacterial_strain=="Db10")) / length(avoid.dat$Worms_ON[which(avoid.dat$Bacterial_strain=="Db10")])
length(which(avoid.dat$Worms_ON+avoid.dat$Worms_OFF+avoid.dat$Excl_from_sub20 <20 & avoid.dat$Bacterial_strain=="OP50")) / length(avoid.dat$Worms_ON[which(avoid.dat$Bacterial_strain=="OP50")])

#create another dataset that only includes N2 when it was the "focal" strain in Block1 (used for predictions)
avoid_onesetN2 = avoid.dat %>% 
  filter(case_when(C.elegans_strain=="N2" ~ as.numeric(Block)<2,
                   TRUE ~ as.numeric(Block)<7))


#### FIGURE 2: C. elegans strains vary in their level of avoidance of the 
#parasite Serratia marcescens (strain Db10). ####
ggplot(data=avoid_onesetN2, 
       aes(x=C.elegans_strain,
           y=Worms_OFF/(Worms_OFF+Worms_ON), fill=Bacterial_strain:Count_time)) + 
  stat_summary(pch=21, size=0.6, fun.data = "mean_se", 
               aes(x=factor(C.elegans_strain, 
                            levels=c("CB4856","LKC34", "JU258", "JU775", "EG4725",
                                     "CX11314","N2","DL238","JT11398", "ED3017","MY16","MY23")),
                            y=Worms_OFF/(Worms_OFF+Worms_ON)), 
               position=position_dodge(width=0.5)) +
  theme_bw() + ylab("Proportion of individuals avoiding bacteria") + xlab("C. elegans strain") +
  scale_fill_manual(values=c("OP50:early" = "lightblue","OP50:late" = "cornflowerblue", "Db10:early" = "pink", "Db10:late" = "coral"),
                    labels=c("OP50-18h","OP50-24h", "Db10-18h","Db10-24h"),
                    name="Bacteria-time") +
  theme(axis.text.x=element_text(angle = 90, hjust=0.95, vjust = 0.5))


###############################################
####### Variation in avoidance behavior #######

##### First: Is the observed behavior pathogen avoidance? #####
# Part A: Do Ce leave lawns of Db10 more than OP50?
# Part B: Is there an interaction between Ce strain and Bact strain to suggest the response is specific to pathogens?

avoid.dat$Plate_ID = paste(avoid.dat$Block, avoid.dat$Plate_number, sep="_")

# Model 1A
#####**************####
model1A = glmer(cbind(Worms_OFF, Worms_ON) ~ C.elegans_strain * Bacterial_strain + Count_time + (1|Plate_ID), 
           data=avoid.dat, family="binomial",
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(model1A) # Table S1
Anova(model1A)

#check for overdispersion with function from Ben Bolker's GLMM FAQ (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(model1A) #ratio = 0.968

#Model 1B
#Model with Block is singular, so I just left plate ID
model1B = glmer(cbind(Worms_OFF, Worms_ON) ~ C.elegans_strain + Bacterial_strain + Count_time + (1|Plate_ID), 
                  data=avoid.dat, family="binomial",
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(model1B) #Table S1
Anova(model1B)

overdisp_fun(model1B) #ratio = 0.73

#compare models
anova(model1A, model1B)
lrtest(model1A, model1B)

#Estimated marginal means from the non-interaction model (in log-odds)
#i.e. what is the mean effect of each Bacterial strain across all dozen Ce strains?
emmeans(model1B, specs = ~Bacterial_strain)

#for OP50 avoidance
plogis(-3.50) #convert log-odds to probability
#CI for emmeans of OP50
plogis(-3.83)
plogis(-3.17)

#Db10 avoidance
plogis(2.21) #convert log-odds to probability
#CI for emmeans of Db10
plogis(1.95)
plogis(2.46)



########### Try more "liberal" interpretation of data - including missing worms in the count as 'off lawn' ##########
#Results are consistent#

avoid.dat$Worms_OFF_20 =(20-avoid.dat$Excl_from_sub20)-avoid.dat$Worms_ON
cbind(avoid.dat$Worms_OFF, avoid.dat$Worms_OFF_20, avoid.dat$Worms_ON)

avoid.dat[which(avoid.dat$Worms_OFF_20<0),]

#Model 2A
model2A = glmer(cbind(Worms_OFF_20, Worms_ON) ~ C.elegans_strain * Bacterial_strain + Count_time + (1|Plate_ID), 
                  data=avoid.dat, family="binomial",
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(model2A) #Table S2
Anova(model2A)

#Model 2B
model2B = glmer(cbind(Worms_OFF_20, Worms_ON) ~ C.elegans_strain + Bacterial_strain + Count_time + (1|Plate_ID), 
                     data=avoid.dat, family="binomial",
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(model2B) #Table S2
Anova(model2B)


emmeans(model2B, specs = ~Bacterial_strain) 
plogis(-2.7)
plogis(2.4)

##### Second: Is there variation between strains just in the level of their responses to Db10? #####
# Is there a main effect of strain?
# Do the strains have different rates of leaving the lawn, i.e. a strain x time interaction? 

####### dataset with just Db10
avoid.Db10 = subset(avoid.dat,Bacterial_strain=="Db10")

#Model 3A - with interaction 
model3A = glmer(data=avoid.Db10, 
                      cbind(Worms_OFF, Worms_ON) ~ C.elegans_strain*Count_time + (1|Block/Plate_ID), 
                      family="binomial")
summary(model3A)
Anova(model3A)
overdisp_fun(model3A) #ratio = 0.58


#Model 3B - without interaction
model3B = glmer(data=avoid.Db10, 
                        cbind(Worms_OFF, Worms_ON) ~ C.elegans_strain + Count_time + (1|Block/Plate_ID), 
                        family="binomial")
summary(model3B)
Anova(model3B)
overdisp_fun(model3B) #ratio = 0.58


#estimated marginal means
emmeans(model3B, specs = ~Count_time)
plogis(1.62)
plogis(2.93)

#compare models; simpler model is better
anova(model3A, model3B)
lrtest(model3A, model3B)

##### Third: What is the broad-sense heritability of avoidance? #####

# 18h time point
model4A = glmer(cbind(Worms_OFF, Worms_ON) ~ 1 + (1|C.elegans_strain), 
                        data=subset(avoid_onesetN2, Bacterial_strain=="Db10" & Count_time=="early"), family="binomial")
summary(model4A)
var(ranef(model4A )[[1]])/ (var(resid(model4A ))+var(ranef(model4A )[[1]])) #0.119

#24h time point
model4B = glmer(cbind(Worms_OFF, Worms_ON) ~ 1 + (1|C.elegans_strain), 
                  data=subset(avoid_onesetN2, Bacterial_strain=="Db10" & Count_time=="late"), family="binomial")
summary(model4B)
var(ranef(model4B )[[1]])/ (var(resid(model4B ))+var(ranef(model4B )[[1]])) #0.255

#combined
model4C = glmer(cbind(Worms_OFF, Worms_ON) ~ 1 + (1|C.elegans_strain), 
                   data=subset(avoid_onesetN2, Bacterial_strain=="Db10"), family="binomial")
summary(model4C)
var(ranef(model4C )[[1]])/ (var(resid(model4C ))+var(ranef(model4C)[[1]])) #0.116

###### For use in later correlation: generate one overall prediction for avoidance,#####
# accounting for both times of counting, with block control
avoid.Db10$avoid_pred = predict(model3B, re.form=~(1|Block), type="response")
avoid.Db10$avoid_pred_LO = predict(model3B, re.form=~(1|Block), type="link")


#get predictions for each strain -- only include N2 in first block for prediction
#also compute mean/se for each strain for comparison
avoid.Db10_summary = avoid.Db10 %>% group_by(Block, C.elegans_strain) %>%
  summarize(avoid_pred= mean(avoid_pred), 
            avoid_pred_LO = mean(avoid_pred_LO), 
            avoid_mean = mean(Worms_OFF/(Worms_ON + Worms_OFF)), 
            se_avoid=(sd(Worms_OFF/(Worms_ON + Worms_OFF))/sqrt(length(Worms_OFF)))) %>%
  filter(case_when(C.elegans_strain=="N2" ~ as.numeric(Block)<2,
                   TRUE ~ as.numeric(Block)<7)) #only consider prediction for N2 in the first block, rather than all the blocks
avoid.Db10_summary


###############################################
####### Variation in fitness #############

#read in the data
popgrowth1 = read.csv("./Data/Population growth assay.csv", header=T)
str(popgrowth1)

#Two plates had very severe burrowing, so they were excluded
popgrowth = popgrowth1 %>% filter(Burrowing_exclude!=1)

#Set variables as factors
popgrowth$Block = factor(popgrowth$Block)
popgrowth$Plate = factor(popgrowth$Plate)

#set N2 as reference strain
popgrowth$Strain = factor(popgrowth$Strain, levels=c("N2", "CX11314", "JU258", "MY23", "CB4856",
                                                     "EG4725", "ED3017", "DL238", "JT11398", 
                                                     "LKC34", "MY16", "JU775"))

#LMM of the mean count of each drop with a fixed effect for Ce strain and 
#random-intercept effect for repeated measures from the same Plate nested within Block

model5 = lmer(Count_20ul ~ Strain + (1|Block/Plate), data=popgrowth, REML=FALSE)
summary(model5)
Anova(model5)


#Fitness model - Generate prediction with random effect of Block included
popgrowth$predict_pgmod1 = predict(model5, re.form=~(1|Block)) 

#Predictions and means at strain level - only use first block of data for N2
pg_summ = popgrowth %>% group_by(Block, Strain) %>%
  summarize(pred_count = mean(predict_pgmod1), mean_count = mean(Count_20ul), 
            se_count=(sd(Count_20ul)/sqrt(length(Count_20ul)))) %>%
  filter(case_when(Strain=="N2" ~ as.numeric(Block)<2,
                   TRUE ~ as.numeric(Block)<6))
pg_summ


######FIGURE 3: C. elegans strains vary in their survival and reproduction in the absence of the parasite. ####
ggplot(data=pg_summ, aes(x=factor(Strain, levels= Strain[order(mean_count)]), y =pred_count)) + 
  geom_errorbar(aes(ymin=mean_count-se_count, ymax=mean_count+se_count), width=0) +
  geom_point(aes(y=mean_count), size=2.5, fill="cornflowerblue", pch=21) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust=0.95, vjust = 0.5)) +
  ylab("Number of individuals in 20ul drop") + xlab("C. elegans strain")


###### Compare fitness estimates to Zhang data #####
zd = read.csv("./Data/Zhang et al 2021 fecundity data.csv", header=T)

#which strains do we have in common
zd %>% filter(strain %in% pg_summ$Strain)

#make a combined dataset to compare
comp = zd %>% inner_join(pg_summ, by=c("strain"="Strain"))
plot(pred_count~mean_b, data=comp)
plot(mean_count~mean_b, data=comp)

summary(lm(pred_count~mean_b, data=comp)) # vs our predictions
summary(lm(mean_count~mean_b, data=comp)) # vs our means

ggplot(data=comp, aes(x=mean_b, y=mean_count)) +
 geom_smooth(method="lm", formula=y~x, col="black") +
  geom_point() + theme_bw() + xlab("Estimate of fitness - Zhang et al. 2021") + 
  ylab("Estimate of fitness - this paper")
  
  
  
###############################################
####### Are avoidance and fitness correlated? #############


##### Combine predictions for avoidance and fitness (with block controls) 
#into one dataset and test for relationship

avoid_fit_comb = avoid.Db10_summary %>% right_join(pg_summ, by=c("C.elegans_strain" = "Strain"))

av.fit.mod = lm(pred_count~avoid_pred_LO, data=avoid_fit_comb)
summary(av.fit.mod) # Table S6
anova(av.fit.mod)


min(avoid_fit_comb$avoid_pred_LO)
max(avoid_fit_comb$avoid_pred_LO)
newdat=data.frame(avoid_pred_LO=seq(1, 3.5, by=0.01))
newdat$avoid_pred=exp(newdat$avoid_pred_LO)/(1+exp(newdat$avoid_pred_LO))
newdat2 = cbind(newdat, predict(av.fit.mod, newdata =  newdat, interval = 'confidence', level=0.95))

##########
#colorblind friendly palette from http://mkweb.bcgsc.ca/colorblind/
  
ggplot(data= newdat2, aes(x=avoid_pred, y=fit)) + 
  geom_line(data=newdat2, aes(x=avoid_pred, y=fit), col="gray40") +
  geom_ribbon(data=newdat2, aes(ymin=lwr,ymax=upr), alpha=0.3, fill="gray70") + 
  geom_errorbar(data= avoid_fit_comb,aes(x=avoid_mean, y= mean_count, xmin=avoid_mean-se_avoid, xmax=avoid_mean+se_avoid, col=C.elegans_strain), width=0, alpha=0.5) +
  geom_errorbar(data= avoid_fit_comb,aes(x=avoid_mean, y= mean_count, ymin=mean_count-se_count, ymax=mean_count+se_count, col=C.elegans_strain), width=0, alpha=0.5) +
  geom_point(data=avoid_fit_comb,aes(x=avoid_mean, y=mean_count, fill=C.elegans_strain), pch=22, size=2, alpha=0.5) +
  geom_point(data= avoid_fit_comb, aes(x=avoid_pred, y= pred_count, fill=C.elegans_strain), pch=21, size=3.5) + 
  scale_fill_manual(values =c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B")) +
  scale_color_manual(values =c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B")) +
  theme_bw() +
  theme(legend.position="none") + 
  xlab("Proporton of individuals avoiding Db10") + 
  ylab("Count in fitness assay") +
  scale_x_continuous(breaks=c(0.75,0.8,0.85,0.9,0.95))
  
# Do we get a relationship between avoidance and Zhang measure of fitness?
avoid_zfit_comp = avoid_fit_comb %>% right_join(comp)

plot(avoid_pred~mean_b, data=avoid_zfit_comp)
summary(lm(mean_b~avoid_pred_LO, data=avoid_zfit_comp))


# is there a relationship bt op50 avoidance and fitness 
avoid.OP50 = subset(avoid.dat, Bacterial_strain=="OP50")
OP50mod = glmer(data=avoid.OP50, 
                        cbind(Worms_OFF, Worms_ON) ~ C.elegans_strain + Count_time + (1|Block/Plate_number), 
                        family="binomial")
summary(OP50mod)
Anova(OP50mod)
#plot(avoid.mod.noint)

avoid.OP50$predict_OP50 = predict(OP50mod, re.form=~(1|Block)) 

avoid.OP50_summary = avoid.OP50 %>% group_by(Block, C.elegans_strain) %>%
  summarize(avoid_OP50_LO= mean(predict_OP50)) %>%
  filter(case_when(C.elegans_strain=="N2" ~ as.numeric(Block)<2,
                   TRUE ~ as.numeric(Block)<7)) #only consider prediction for N2 in the first block, rather than all the blocks
avoid.OP50_summary
OP50_fit_comb = avoid.OP50_summary %>% right_join(avoid_fit_comb)

summary(lm(pred_count~avoid_OP50_LO, data=OP50_fit_comb)) #table S7
summary(lm(avoid_pred_LO~avoid_OP50_LO, data=OP50_fit_comb)) #table S8


# Citations for R, RStudio, packages 
citation()
RStudio.Version()
citation("tidyverse")
citation("emmeans")
citation("lme4")
citation("car")
citation("lmtest")
citation("stats")
citation("lmerTest")


#>  sessionInfo()
#R version 4.3.0 (2023-04-21)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS Monterey 12.6.7

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: America/New_York
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] pbkrtest_0.5.2  lmerTest_3.1-3  lmtest_0.9-40   zoo_1.8-12      car_3.1-2      
#[6] carData_3.0-5   lme4_1.1-33     Matrix_1.5-4    emmeans_1.8.7   lubridate_1.9.2
#[11] forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2     purrr_1.0.1     readr_2.1.4    
#[16] tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.3        httr2_0.2.3         gh_1.4.0            lattice_0.21-8     
#[5] numDeriv_2016.8-1.1 tzdb_0.4.0          vctrs_0.6.2         tools_4.3.0        
#[9] generics_0.1.3      parallel_4.3.0      curl_5.0.0          fansi_1.0.4        
#[13] pkgconfig_2.0.3     lifecycle_1.0.3     farver_2.1.1        compiler_4.3.0     
#[17] credentials_1.3.2   munsell_0.5.0       sys_3.4.1           usethis_2.2.2      
#[21] pillar_1.9.0        nloptr_2.0.3        crayon_1.5.2        MASS_7.3-58.4      
#[25] openssl_2.0.6       boot_1.3-28.1       abind_1.4-5         nlme_3.1-162       
#[29] tidyselect_1.2.0    mvtnorm_1.2-2       stringi_1.7.12      labeling_0.4.2     
#[33] splines_4.3.0       rprojroot_2.0.3     grid_4.3.0          colorspace_2.1-0   
#[37] cli_3.6.1           magrittr_2.0.3      utf8_1.2.3          broom_1.0.4        
#[41] withr_2.5.0         backports_1.4.1     scales_1.2.1        rappdirs_0.3.3     
#[45] estimability_1.4.1  timechange_0.2.0    gitcreds_0.1.2      askpass_1.1        
#[49] hms_1.1.3           coda_0.19-4         mgcv_1.8-42         rlang_1.1.1        
#[53] gert_1.9.3          Rcpp_1.0.10         xtable_1.8-4        glue_1.6.2         
#[57] rstudioapi_0.14     minqa_1.2.5         jsonlite_1.8.4      R6_2.5.1           
#[61] fs_1.6.2           
 

