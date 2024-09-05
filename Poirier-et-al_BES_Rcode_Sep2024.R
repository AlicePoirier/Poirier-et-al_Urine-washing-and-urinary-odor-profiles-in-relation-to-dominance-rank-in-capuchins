#
##
### Supporting R code for publication in Behavioral Ecology and Sociobiology
### "Dominance effects on urine washing and urinary odor profiles 
### in wild male capuchin monkeys (Cebus imitator)"
### Authors: Alice C. Poirier, Nelle K. Kulick, Suheidy Romero Morales, 
### Marlen Kücklich, Brigitte M. Weiß, Claudia Birkemeyer, Anja Widdig, 
### Amanda D. Melin, and Katharine M. Jack
##
#


# Libraries (here using R 4.4.1)
library(devtools)
library(reshape2)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(vegan)
library(DHARMa)
library(performance)
library(predictmeans)
library(randomForest)




### 1. Prediction 1: UW Behaviour ----
## 1.1. Dataset ----
# Behaviour dataset: Poirier-et-al_BES_Dataset_FocalBehaviour.csv
data.beh<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.beh[,c(1:10)]<- lapply(data.beh[,c(1:10)], factor)
str(data.beh) # 2049 focalnames


range(data.beh$AveRain[which(data.beh$WetDry=="DRY")])
range(data.beh$AveRain[which(data.beh$WetDry=="WET")])


## 1.2. Summary UW per individual, and boxplot UW by rank and season  ----
dev.new()

data.beh.indiv.count<- data.beh %>% group_by(WetDry, Rank, Indiv, Behav) %>% tally()
data.beh.indiv.count<- data.beh.indiv.count[!is.na(data.beh.indiv.count$Behav),]
data.beh.indiv.count[,c(1:4)]<- lapply(data.beh.indiv.count[,c(1:4)], factor)
data.beh.indiv.count


# Figure 2
ggplot(data.beh.indiv.count, aes(x=WetDry, y=n, fill=Rank))+
  geom_boxplot()+
  scale_fill_grey()+
  xlab("Season") +
  ylab("Number of UW events")+
  guides(fill=F)+
  theme_classic(base_size=12)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))+
  theme(legend.justification=c(0,1), legend.position="top")



## 1.3. Behaviour models ----
# 1.3.1. UW Frequency dataset ----
data.beh.count<- data.frame(data.beh %>% 
                              group_by(Beh_type, WetDry, Date, Year, Beh_focalname, Rank, Group, Indiv, z.Age, z.AveRain, z.MeanTemp, z.DayRain) %>% 
                              reframe(count=n()))
data.beh.count$count[which(is.na(data.beh.count$Beh_type))] <- 0
str(data.beh.count)


# 1.3.2. Full and Null behaviour models ----
beh.Full.mod<- glmmTMB(count~Rank*z.AveRain+Group+z.Age+z.MeanTemp*z.DayRain+(1|Indiv)+(1|Date), 
                       data=data.beh.count, ziformula=~1, family=poisson(link="log"))
summary(beh.Full.mod)


beh.Null.mod<- glmmTMB(count~(1|Indiv)+(1|Date), 
                       data=data.beh.count, ziformula=~1, family=poisson(link="log"))
summary(beh.Null.mod)


# 1.3.3. Full-Null model comparison ----
anova(beh.Null.mod, beh.Full.mod, test="Chisq") #***




### 2. Prediction 2: Urine chemical profiles ----
## 2.1. Datasets ----
# 2.1.1 Read csv ----
# Sample List: Poirier-et-al_BESr_Dataset_SampleList.csv
data.spl<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.spl[,c(1:12)]<- lapply(data.spl[,c(1:12)], factor)
str(data.spl) # 147 files

# Compound List: Poirier-et-al_BES_Dataset_CompoundList.csv
data.chem<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
data.chem[,c(1:13)]<- lapply(data.chem[,c(1:13)], factor)
str(data.chem) # 147 files, 525 peaks

# mean + sd number of compounds per sample
Indiv_comp<- as.data.frame(xtabs(~Filename, data.chem))
mean(Indiv_comp$Freq)
sd(Indiv_comp$Freq)



# 2.1.2. Remove rows with NAs as fixed predictors ----
# Sample List
data.spl2<- data.spl[!is.na(data.spl$Rank),]
data.spl2<- data.spl2[!is.na(data.spl2$UW),]
data.spl2<- data.spl2[!is.na(data.spl2$Catch),]
data.spl2<- data.spl2[!is.na(data.spl2$Group),]
data.spl2<- data.spl2[!is.na(data.spl2$meanSPG),]
data.spl2<- data.spl2[!is.na(data.spl2$AveRain),]
data.spl2<- data.spl2[!is.na(data.spl2$DayRain),]
data.spl2<- data.spl2[!is.na(data.spl2$Uvol),]
data.spl2[,c(1:12)]<- lapply(data.spl2[,c(1:12)], factor)
str(data.spl2) # 138 filenames

# Compound List
data.chem2<- data.chem[!is.na(data.chem$Rank),]
data.chem2<- data.chem2[!is.na(data.chem2$UW),]
data.chem2<- data.chem2[!is.na(data.chem2$Catch),]
data.chem2<- data.chem2[!is.na(data.chem2$Group),]
data.chem2<- data.chem2[!is.na(data.chem2$meanSPG),]
data.chem2<- data.chem2[!is.na(data.chem2$AveRain),]
data.chem2<- data.chem2[!is.na(data.chem2$DayRain),]
data.chem2<- data.chem2[!is.na(data.chem2$Uvol),]
data.chem2[,c(1:13)]<- lapply(data.chem2[,c(1:13)], factor)
str(data.chem2) # 14711 obs, 138 filenames, 525 peaks



# 2.1.3. RPA standardization ----
head(data.chem2)
data.chem2<- data.chem2 %>% arrange(Filename, Area)
data.chem.RPA<- data.frame(data.chem2 %>%
                             group_by(Filename) %>% 
                             reframe(Area=Area, SArea=sum(Area), meanArea=sum(Area)/Ncomp, 
                                     RPA=round((Area/sum(Area))*100, digits=4),
                                     logRPA=round(log(((Area/sum(Area))*100)+1), digits=4), 
                                     asinRPA=round((log(asin(sqrt(((Area/sum(Area))*100)/100))+0.01)), digits=4)))
data.chem.RPA<- data.chem.RPA %>% arrange(Filename, Area)

data.chem.std<- cbind(data.chem2, data.chem.RPA[,-(1:2)])
str(data.chem.std)


# 2.1.4. Matrix of samples x compounds ----
data.chem.mat<- xtabs(logRPA~Filename+Compound, data.chem.std)
data.chem.vdim<-vegdist(data.chem.mat, method="bray")



## 2.2. PerMANOVA approach ----
# 2.2.1. Assumption = homogeneity of multivariate dispersion ----
bet.Rank<- betadisper(data.chem.vdim, data.spl2$Rank)
permutest(bet.Rank, permutations=999) #NS

bet.Ind<- betadisper(data.chem.vdim, data.spl2$Indiv)
permutest(bet.Ind, permutations=999) #NS

bet.Grp<- betadisper(data.chem.vdim, data.spl2$Group)
permutest(bet.Grp, permutations=999) #NS

bet.re<- betadisper(data.chem.vdim, data.spl2$Res)
permutest(bet.re, permutations=999) #***

bet.Batch<- betadisper(data.chem.vdim, data.spl2$Batch)
permutest(bet.Batch, permutations=999) #ns

bet.UW<- betadisper(data.chem.vdim, data.spl2$UW)
permutest(bet.UW, permutations=999) #NS

bet.Catch<- betadisper(data.chem.vdim, data.spl2$Catch)
permutest(bet.Catch, permutations=999) #NS


# 2.2.2. Adonis ----
data.spl2<-arrange(data.spl2, Filename)

# Using Indiv as stratum, full-null comparison
Adonis.mod<- adonis2(data.chem.vdim~Rank+z.AveRain+z.Age+Group+UW+Catch+z.SPG+z.Uvol+z.Temp*z.DayRain+Res+Batch, 
                     strata=data.spl2$Indiv, data=data.spl2, permutations=999, method="bray", by=NULL)
Adonis.mod


## Full model using margin to see marginal effects of each predictor
Adonis.mod.marg<- adonis2(data.chem.vdim~Rank+z.AveRain+z.Age+Group+UW+Catch+z.SPG+z.Uvol+z.Temp*z.DayRain+Res+Batch, 
                          strata=data.spl2$Indiv, data=data.spl2, permutations=999, method="bray", by='margin')
Adonis.mod.marg



## 2.3. Random Slope Models approach ----
# 2.3.1. Check that predictors are not too unbalanced ----
# (= one type entirely missing or with only 1-2 occurrences)
xtabs(~Batch+Rank, data.chem.std)    # OK
xtabs(~Batch+Group, data.chem.std)   # several 0 but OK
xtabs(~Batch+UW, data.chem.std)      # OK
xtabs(~Batch+Catch, data.chem.std)   # several 0 but OK

xtabs(~Indiv+UW, data.chem.std)      # several 0 but OK
xtabs(~Indiv+Catch, data.chem.std)   # many 0 -> remove slope

xtabs(~Res+Rank, data.chem.std)      # OK
xtabs(~Res+Group, data.chem.std)     # many 0 -> remove slope
xtabs(~Res+UW, data.chem.std)        # OK
xtabs(~Res+Catch, data.chem.std)     # many 0 -> remove slope


# 2.3.2. Full and Null models ----
head(data.chem.std)

# Null Model
Null.mod=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename) +(1|Compound)
              +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
              data=data.chem.std, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

summary(Null.mod)

# Full Model - long to run !!
Full.mod=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
              +(1+Rank*z.AveRain+z.Age+z.SPG+z.Uvol+z.Temp*z.DayRain||Compound)
              +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
              +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
              data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

summary(Full.mod)


# 2.3.3. Reduced models ----
# long to run !!
# omitting Rank slope
Red.mod.rk=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Rank*seasonality slope
Red.mod.rs=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Age slope
Red.mod.ag=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Group slope
Red.mod.gp=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting AveRain slope
Red.mod.ar=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting UW slope
Red.mod.uw=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Catch slope
Red.mod.ca=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) 
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting SPG slope
Red.mod.sp=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.Uvol+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Uvol slope
Red.mod.vo=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+Rank*z.AveRain+z.Temp*z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Temp slope
Red.mod.te=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting DayRain slope
Red.mod.dr=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))

# omitting Temp*DayRain slope
Red.mod.td=lmer(asinRPA ~ Rank*z.AveRain +z.Age +Group +UW +Catch +z.SPG +z.Uvol +z.Temp*z.DayRain +(1|Batch) +(1|Indiv) +(1|Res) +(1|Filename)
                +(1+z.Age+z.SPG+z.Uvol+Rank*z.AveRain+z.Temp+z.DayRain||Compound)
                +(1|Compound:Rank) +(1|Compound:Group) +(1|Compound:UW) +(1|Compound:Catch)
                +(1|Batch:Rank) +(1|Batch:Group) +(1|Batch:UW) +(1|Batch:Catch) +(1|Indiv:UW) +(1|Res:Rank) +(1|Res:UW),
                data=data.chem.std, verbose=10, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))


# 2.3.4. Anovas
anova(Null.mod, Full.mod, test="Chisq") #

anova(Red.mod.rk, Full.mod, test="Chisq") #
anova(Red.mod.rs, Full.mod, test="Chisq") #
anova(Red.mod.ag, Full.mod, test="Chisq") #
anova(Red.mod.gp, Full.mod, test="Chisq") #
anova(Red.mod.ar, Full.mod, test="Chisq") #
anova(Red.mod.uw, Full.mod, test="Chisq") #
anova(Red.mod.ca, Full.mod, test="Chisq") #
anova(Red.mod.sp, Full.mod, test="Chisq") #
anova(Red.mod.vo, Full.mod, test="Chisq") #
anova(Red.mod.te, Full.mod, test="Chisq") #
anova(Red.mod.dr, Full.mod, test="Chisq") #
anova(Red.mod.td, Full.mod, test="Chisq") #





### 3. Check linear model assumptions ----
# Choose correct model name
#mymod<- beh.Full.mod  
mymod<- Full.mod  
mymod<- beh.Full.mod  

## 3.1. Package predictmeans ----
residplot(mymod)

## 3.2. Package performance ----
check_collinearity(mymod)
check_distribution(mymod)
check_heteroscedasticity(mymod)
check_autocorrelation(mymod)
#check_homogeneity(mymod) #bugged
check_singularity(mymod)
check_outliers(mymod)



# end of script

