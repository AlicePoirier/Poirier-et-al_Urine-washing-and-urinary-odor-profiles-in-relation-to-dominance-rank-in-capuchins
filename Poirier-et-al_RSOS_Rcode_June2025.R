#
##
### Supporting R code for publication in Royal Society Open Science:
###
### "Urine washing and urinary odor profiles in relation to dominance  
### rank status in wild male capuchin monkeys (Cebus imitator)"
###
### Authors: Alice C. Poirier, Nelle K. Kulick, Suheidy Romero Morales, 
### Marlen Kücklich, Brigitte M. Weiß, Claudia Birkemeyer, Anja Widdig, 
### Amanda D. Melin, and Katharine M. Jack
##
#

# Libraries (here using R 4.5.0)
library(devtools)
library(reshape2)
library(lubridate)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(sjPlot)
library(vegan)
library(DHARMa)
library(predictmeans)



############
### 1. Prediction 1: UW Behaviour ----
## 1.1. Datasets ----
# 1.1.1. Behaviour dataset: Poirier-et-al_RSOS_Dataset_FocalBehaviour_June2025.csv ----
data.beh<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.beh[,c(1:13)]<- lapply(data.beh[,c(1:13)], factor)
str(data.beh) # 2641 focal unique IDs


# 1.1.2. UW Frequency dataset ----
# Count 1 for UW behav
data.beh2<- mutate(data.beh, UW_count=rep(0))
data.beh2$UW_count[which(data.beh2$UW_behav=="UW")] <- 1
head(data.beh2)

data.beh.count<- data.frame(data.beh2 %>% 
                              group_by(Focal_uniqueID, Date, Male, Rank, Group, Observer, Focal_duration, Age, Seasonality) %>% 
                              reframe(UW_count=sum(UW_count)))
data.beh.count[,c(1:7)]<- lapply(data.beh.count[,c(1:7)], factor)
head(data.beh.count)


# Transform continuous predictors
data.beh.count$z.Age<- as.vector(scale(data.beh.count$Age))
data.beh.count$z.Seasonality<- as.vector(scale(data.beh.count$Seasonality))

# Transform focal duration to seconds
data.beh.count$Duration_sec <- as.numeric(hms(data.beh.count$Focal_duration))
head(data.beh.count)



## 1.2. Descriptive stats on focal observations ----
# 1.2.1. Age range ----
range(data.beh.count$Age)
mean(data.beh.count$Age)
sd(data.beh.count$Age)


# 1.2.2. Total observation time ----
(sum(data.beh.count$Duration_sec))/3600 #hours

(range(data.beh.count$Duration_sec))/60
(mean(data.beh.count$Duration_sec))/60
(sd(data.beh.count$Duration_sec))/60


# 1.2.3. Frequency of UW per focal ----
sum(data.beh.count$UW_count)
mean(data.beh.count$UW_count)
sd(data.beh.count$UW_count)


# 1.2.4. Descriptive stats per male ----
data.beh.count.male<- data.frame(data.beh.count %>% group_by(Male) %>% 
                                     reframe(Focal_duration_hour=round(sum(Duration_sec)/3600, 1),
                                             UW_count_total=sum(UW_count)))
data.beh.count.male<- mutate(data.beh.count.male, UW_hour_rate=round(UW_count_total/Focal_duration_hour,2))
data.beh.count.male

range(data.beh.count.male$Focal_duration_hour)
mean(data.beh.count.male$Focal_duration_hour)
sd(data.beh.count.male$Focal_duration_hour)

range(data.beh.count.male$UW_hour_rate)
mean(data.beh.count.male$UW_hour_rate)
sd(data.beh.count.male$UW_hour_rate)


# 1.2.5. Descriptive stats per rank status and male ----
data.beh.count.rank<- data.frame(data.beh.count %>% group_by(Rank, Male) %>% 
                                   reframe(Focal_duration_hour=round(sum(Duration_sec)/3600, 1),
                                           UW_count_total=sum(UW_count)))
data.beh.count.rank<- mutate(data.beh.count.rank, UW_hour_rate=round(UW_count_total/Focal_duration_hour,2))
data.beh.count.rank

range(data.beh.count.rank[data.beh.count.rank$Rank=="alpha",]$UW_hour_rate)
mean(data.beh.count.rank[data.beh.count.rank$Rank=="alpha",]$UW_hour_rate)
sd(data.beh.count.rank[data.beh.count.rank$Rank=="alpha",]$UW_hour_rate)

range(data.beh.count.rank[data.beh.count.rank$Rank=="sub",]$UW_hour_rate)
mean(data.beh.count.rank[data.beh.count.rank$Rank=="sub",]$UW_hour_rate)
sd(data.beh.count.rank[data.beh.count.rank$Rank=="sub",]$UW_hour_rate)


# 1.2.6. Descriptive stats per group ID, rank status and male (Suppl. Table S1) ----
data.beh.count.group<- data.frame(data.beh.count %>% group_by(Group, Rank, Male) %>% 
                                   reframe(Focal_duration_hour=round(sum(Duration_sec)/3600, 1),
                                           UW_count_total=sum(UW_count)))
data.beh.count.group<- mutate(data.beh.count.group, UW_hour_rate=round(UW_count_total/Focal_duration_hour,2))
data.beh.count.group




## 1.3. Behaviour models ----
# 1.3.1. Full and Null behaviour models ----
head(data.beh.count)
data.beh.count2<- na.omit(data.beh.count)

# Full model including the interaction Rank*z.Seasonality
beh.Full.mod.int<- glmmTMB(UW_count ~ Rank*z.Seasonality +z.Age +Group +(1|Male) +(1|Date) +(1|Observer) +offset(log(Duration_sec)), 
                       data=data.beh.count2, ziformula=~1, family=poisson(link="log"))
summary(beh.Full.mod.int)

# Reduced model without the interaction Rank*z.Seasonality
beh.Full.mod<- glmmTMB(UW_count ~ Rank +z.Seasonality +z.Age +Group +(1|Male) +(1|Date) +(1|Observer) +offset(log(Duration_sec)), 
                       data=data.beh.count2, ziformula=~1, family=poisson(link="log"))
summary(beh.Full.mod)

# Null model
beh.Null.mod<- glmmTMB(UW_count ~ 1 +(1|Male) +(1|Date) +(1|Observer) +offset(log(Duration_sec)), 
                       data=data.beh.count2, ziformula=~1, family=poisson(link="log"))
summary(beh.Null.mod)



# 1.3.2. Model comparisons ----
anova(beh.Full.mod, beh.Full.mod.int, test="Chisq") 
#* significant so important to keep the interaction term

anova(beh.Null.mod, beh.Full.mod.int, test="Chisq") 
#***

# Full model estimates
Anova(beh.Full.mod.int, type="II", error.estimate="pearson")




## 1.4. Figure 2 ----
# Plot model
set_theme(base = theme_classic(base_size=12))
plot_model(beh.Full.mod.int, type="pred", terms=c("z.Seasonality", "Rank"), condition=c(Duration_sec=1),
           axis.title=c("Seasonality (z-transformed)", "Predicted UW count values"), bias_correction=T,
           title="", show.legend=F, line.size=1, show.data=T, jitter=0.0003)



####
### 2. Prediction 2: Urine chemical profiles ----
## 2.1. Datasets ----
# 2.1.1 Sample list ----
# Sample List: Poirier-et-al_RSOS_Dataset_GCMS_SampleList_June2025.csv
data.spl<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.spl[,c(1:9)]<- lapply(data.spl[,c(1:9)], factor)
str(data.spl) # 153 files

# z-transform variables
data.spl$z.Age<- as.vector(scale(data.spl$Age))
data.spl$z.Volume<- as.vector(scale(data.spl$Volume))
data.spl$z.SPG<- as.vector(scale(data.spl$SPG))
data.spl$z.Temp<- as.vector(scale(data.spl$Temp))
data.spl$z.Seasonality<- as.vector(scale(data.spl$Seasonality))

# Add continuous time
data.spl<-arrange(data.spl, Date)
data.spl<- mutate(data.spl, n.Date=as.numeric(Date))

# Number of samples per male, group and rank (Suppl Table S1)
data.spl.maleGR<- data.frame(data.spl %>% group_by(Male, Group, Rank) %>% tally())
data.spl.maleGR

# Number of samples per male
data.spl.male<- data.frame(data.spl %>% group_by(Male) %>% tally())
mean(data.spl.male$n)
sd(data.spl.male$n)

# Urine volume
range(data.spl$Volume)
mean(data.spl$Volume)
sd(data.spl$Volume)



# 2.1.2. Compound list ----
# Compound List: Poirier-et-al_RSOS_Dataset_GCMS_CompoundList_June2025.csv
data.chem<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.chem[,c(1:2)]<- lapply(data.chem[,c(1:2)], factor)
str(data.chem) 

# Number of compounds per sample
Indiv_comp<- as.data.frame(xtabs(~Filename, data.chem))
names(Indiv_comp)[names(Indiv_comp)=="Freq"] <- "Ncomp"
mean(Indiv_comp$Ncomp)
sd(Indiv_comp$Ncomp)

data.chem2 <- merge(data.chem, Indiv_comp, by="Filename")
data.chem2[,c(1:2)]<- lapply(data.chem2[,c(1:2)], factor)
str(data.chem2)



# 2.1.3. Merge compound list and sample list ----
data.spl.chem <- merge(data.spl, data.chem2, by="Filename")
data.spl.chem[,c(1:9,21)]<- lapply(data.spl.chem[,c(1:9,21)], factor)
str(data.spl.chem) # 153 files, 88 compounds



# 2.1.4. RPA standardization ----
data.spl.chem<- data.spl.chem %>% arrange(Filename, Area)
data.chem.RPA<- data.frame(data.spl.chem %>%
                             group_by(Filename) %>% 
                             reframe(Area=Area, SArea=sum(Area), meanArea=sum(Area)/Ncomp, 
                                     RPA=round((Area/sum(Area))*100, digits=4),
                                     logRPA=round(log(((Area/sum(Area))*100)+1), digits=4), 
                                     asinRPA=round((log(asin(sqrt(((Area/sum(Area))*100)/100))+0.01)), digits=4)))
data.chem.RPA<- data.chem.RPA %>% arrange(Filename, Area)

data.chem.std<- cbind(data.spl.chem, data.chem.RPA[,-(1:2)])
str(data.chem.std)


# 2.1.5. Matrix of samples x compounds ----
data.chem.mat<- xtabs(logRPA~Filename+Comp, data.chem.std)
data.chem.vdim<-vegdist(data.chem.mat, method="bray")

data.chem.nmds<- metaMDS(data.chem.vdim, k=3, autotransform=F, noshare=F, wascores=F, expand=T, trymax=200, trace=F)
data.chem.nmds.sco<- vegan::scores(data.chem.nmds)
data.chem.nmds$stress    #0.15
meanNMDS<- rowMeans(data.chem.nmds.sco)

data.spl<- data.spl %>% arrange(Filename)
data.chem.res<- cbind(data.spl, data.chem.nmds.sco, meanNMDS)
data.chem.res<-arrange(data.chem.res, Filename)
head(data.chem.res)



## 2.2. PerMANOVA approach ----
# 2.2.1. Assumption = homogeneity of multivariate dispersion ----
bet.Rank<- betadisper(data.chem.vdim, data.chem.res$Rank)
permutest(bet.Rank, permutations=999) #NS

bet.Grp<- betadisper(data.chem.vdim, data.chem.res$Group)
permutest(bet.Grp, permutations=999) #**

bet.Catch<- betadisper(data.chem.vdim, data.chem.res$Catch)
permutest(bet.Catch, permutations=999) #NS

bet.Mal<- betadisper(data.chem.vdim, data.chem.res$Male)
permutest(bet.Mal, permutations=999) #***

bet.Han<- betadisper(data.chem.vdim, data.chem.res$Handler)
permutest(bet.Han, permutations=999) #***

bet.Batch<- betadisper(data.chem.vdim, data.chem.res$Batch)
permutest(bet.Batch, permutations=999) #***



# 2.2.2. Adonis ----
data.chem.res<-arrange(data.chem.res, Filename)

## Using Indiv as stratum, full-null comparison
# without interaction term rank*seasonality
Adonis.mod<- adonis2(data.chem.vdim ~ Rank+z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +Handler +Batch +n.Date, 
                     strata=data.chem.res$Male, data=data.chem.res, permutations=999, method="bray", by=NULL)
Adonis.mod #***


## Full model using margin to test the marginal effect of each predictor
# with interaction term rank*seasonality
Adonis.mod.marg.int<- adonis2(data.chem.vdim ~ Rank*z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +Handler +Batch +n.Date, 
                          strata=data.chem.res$Male, data=data.chem.res, permutations=999, method="bray", by='margin')
Adonis.mod.marg.int 
#the interaction term is non-significant so we can use model without interaction


# without interaction term rank*seasonality
Adonis.mod.marg<- adonis2(data.chem.vdim ~ Rank+z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +Handler +Batch +n.Date, 
                          strata=data.chem.res$Male, data=data.chem.res, permutations=999, method="bray", by='margin')
Adonis.mod.marg




## 2.3. Random Slope Models approach ----
# 2.3.1. Check that predictors are not too unbalanced ----
# (= one type entirely missing or with only 1-2 occurrences)
xtabs(~Batch+Rank, data.chem.std)       # OK
xtabs(~Batch+Group, data.chem.std)      # OK
xtabs(~Batch+Catch, data.chem.std)      # OK

xtabs(~Male+Catch, data.chem.std)       # ok

xtabs(~Handler+Rank, data.chem.std)     # OK
xtabs(~Handler+Group, data.chem.std)    # many 0 -> remove slope
xtabs(~Handler+Catch, data.chem.std)    # OK

# Check VIF
mod.vif=lm(asinRPA ~ Rank+z.Seasonality+Catch+Group+z.Age+z.SPG+z.Volume+z.Temp+n.Date, data=data.chem.std)
vif(mod.vif) # ok, max vif = 2.4 for group



# 2.3.2. Full and Null models ----
head(data.chem.std)

## Null models
# with interaction term rank*seasonality
Null.mod.int=lmer(asinRPA ~ Rank*z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +n.Date +(1|Male) +(1|Date) +(1|Handler) +(1|Batch) +(1|Filename)
                  +(1+z.SPG+z.Volume+z.Temp+n.Date||Comp) +(1|Comp:Catch)
                  +(1|Male:Catch) +(1|Handler:Rank) +(1|Handler:Catch) +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:Catch),
                  data=data.chem.std, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Null.mod.int)

# without interaction term rank*seasonality
Null.mod=lmer(asinRPA ~ Rank+z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +n.Date +(1|Male) +(1|Date) +(1|Handler) +(1|Batch) +(1|Filename)
              +(1+z.SPG+z.Volume+z.Temp+n.Date||Comp) +(1|Comp:Catch)
              +(1|Male:Catch) +(1|Handler:Rank) +(1|Handler:Catch) +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:Catch),
              data=data.chem.std, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Null.mod)

anova(Null.mod, Null.mod.int, test="Chisq") #NS
# non-significant test so we can omit the interaction term


# Full model without interaction term
Full.mod=lmer(asinRPA ~ Rank +z.Seasonality +z.Age +Group +Catch +z.SPG +z.Volume +z.Temp +n.Date +(1|Male) +(1|Date) +(1|Handler) +(1|Batch) +(1|Filename)
                  +(1+z.Seasonality+z.Age+z.SPG+z.Volume+z.Temp+n.Date||Comp)
                  +(1|Comp:Rank) +(1|Comp:Group) +(1|Comp:Catch)
                  +(1|Male:Catch) +(1|Handler:Rank) +(1|Handler:Catch) +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:Catch),
                  data=data.chem.std, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.mod)

anova(Null.mod, Full.mod, test="Chisq") #NS





### 3. Check linear model assumptions ----
# Choose correct model name
#mymod<- beh.Full.mod  
mymod<- Full.mod


## 3.1. Package predictmeans ----
residplot(mymod)

## 3.2. Package DHARMa ----
simulationOutput<- simulateResiduals(fittedModel=mymod, n=1000, refit=F, use.u=F)
testResiduals(simulationOutput) 




#save.image(file="current_RSOS.RData")

# end of script

