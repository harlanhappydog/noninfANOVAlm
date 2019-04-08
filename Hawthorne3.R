library(gee)
library(geepack)
library(BayesFactor)
library(RCurl)
library(httr)

#############################################################
# This loads custom code for ANOVA non-inferiority testing from github repo:
script <- getURL("https://raw.githubusercontent.com/harlanhappydog/noninfANOVAlm/master/noninfANOVA.R", ssl.verifypeer = FALSE)
eval(parse(text = script))


#############################################################
## This reads in the data, and formats it appropriately:
Hdata <-read.csv(text= getURL("https://raw.githubusercontent.com/harlanhappydog/noninfANOVAlm/master/analysisdataset2018.csv"))

Hdata$group<-as.factor(Hdata$group)
Hdata$t<-as.factor(Hdata$t)
Hdata$participant_ID <-as.factor(Hdata$participant_ID)

Hdata<- Hdata[!(Hdata$group)=="",]
Hdata$group<-as.factor(as.character(Hdata$group))

Hdata<-Hdata[order(Hdata$participant_ID),]


#############################################################
## This reads reduces the data to subset of complete cases:

Hdata<-Hdata[Hdata$itt==1,]

# option to select a random subset of 1000 samples from the dataset:
# Hdata <-Hdata[sample(1:dim(Hdata)[1], 1000),]

Hdata_complete<-Hdata[Hdata$participant_ID%in%(c(names(table(Hdata$participant_ID)[table(Hdata$participant_ID)==2]))),]

Hdata_complete$participant_ID<-as.factor(as.character((Hdata_complete$participant_ID)))
Hdata_complete<-Hdata_complete[order(Hdata_complete$participant_ID),]


#############################################################
## This creates a "wide" version of the complete cases data:

base_data<-Hdata_complete[Hdata_complete$t=="Baseline",][,c("participant_ID","DrinkingDays", "drinksperday", "totaldrinking", "group")]

followup_data<-Hdata_complete[Hdata_complete$t=="Followup",][,c("participant_ID","DrinkingDays", "drinksperday", "totaldrinking")]

side_data<-merge(base_data , followup_data, by="participant_ID")
side_data$DrinkingDays.diff<-side_data$DrinkingDays.y-side_data$DrinkingDays.x
side_data$drinksperday.diff<-side_data$drinksperday.y-side_data$drinksperday.x
side_data$totaldrinking.diff<-side_data$totaldrinking.y-side_data$totaldrinking.x


#############################################################
## This analysis code reproduces the results in the paper:

## DrinkingDays
Hdata$group<-relevel(Hdata$group,"A")
summary(geeglm(DrinkingDays ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))
Hdata$group<-relevel(Hdata$group,"C")
summary(geeglm(DrinkingDays ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))


## totaldrinking
Hdata$group<-relevel(Hdata$group,"A")
summary(geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))
Hdata$group<-relevel(Hdata$group,"C")
summary(geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))


## drinksperday
Hdata$group<-relevel(Hdata$group,"A")
summary(geeglm(drinksperday ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))
Hdata$group<-relevel(Hdata$group,"C")
summary(geeglm(drinksperday ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata))

# Note that there is a typo in paper:  0.86 should be 0.89 and vice-versa

#############################################################
# Here are three alternative ways to do the analysis (DrinkingDays):
side_data$group<-relevel(side_data $group,"A")
Hdata$group<-relevel(Hdata$group,"A")

summary(geeglm(DrinkingDays ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(DrinkingDays.y  ~ group + DrinkingDays.x , data= side_data))
summary(lm(DrinkingDays.diff  ~ group , data= side_data))


side_data$group<-relevel(side_data $group,"C")
Hdata$group<-relevel(Hdata$group,"C")

summary(geeglm(DrinkingDays ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(DrinkingDays.y  ~ group + DrinkingDays.x , data= side_data))
summary(lm(DrinkingDays.diff  ~ group , data= side_data))


#############################################################
# Here are three alternative ways to do the analysis (totaldrinking):
side_data$group<-relevel(side_data $group,"A")
Hdata$group<-relevel(Hdata$group,"A")

summary(geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(totaldrinking.y  ~ group + totaldrinking.x , data= side_data))
summary(lm(totaldrinking.diff  ~ group , data= side_data))


side_data$group<-relevel(side_data $group,"C")
Hdata$group<-relevel(Hdata$group,"C")

summary(geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(totaldrinking.y  ~ group + totaldrinking.x , data= side_data))
summary(lm(totaldrinking.diff  ~ group , data= side_data))

#############################################################
# Here are three alternative ways to do the analysis (drinksperday):
side_data$group<-relevel(side_data $group,"A")
Hdata$group<-relevel(Hdata$group,"A")

summary(geeglm(drinksperday ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(drinksperday.y  ~ group + drinksperday.x , data= side_data))
summary(lm(drinksperday.diff  ~ group , data= side_data))


side_data$group<-relevel(side_data $group,"C")
Hdata$group<-relevel(Hdata$group,"C")

summary(geeglm(drinksperday ~ group*t + group+t, id= participant_ID, corstr="exchangeable", data= Hdata))
summary(lm(drinksperday.y  ~ group + drinksperday.x , data= side_data))
summary(lm(drinksperday.diff  ~ group , data= side_data))


#############################################################
# Here is how we recommend doing the non-inferiority analysis
# in a linear regression context (drinksperday):

Xmatrix <- model.matrix(drinksperday.diff  ~ group, data= side_data)
lmmodel <- lm(drinksperday.diff  ~ group , data= side_data)

R2 <- summary(lmmodel)$r.squared
Fstat <- summary(lmmodel)$fstatistic[1]
K <- dim(Xmatrix)[2] - 1
N <- dim(Xmatrix)[1]
Delta <- 0.01
 
pf(Fstat,df1=K,df2=N-K-1,ncp=(N*Delta)/(1-Delta),lower.tail=TRUE)

# Here is an analogous Bayesian testing scheme to consider:
linearReg.R2stat(N=N, p=K, R2= R2, simple=TRUE)

#############################################################
# Here is the equivalent non-inferiority analysis
# in an ANOVA (equal variances) (drinksperday):

Xmatrix <- model.matrix(drinksperday.diff  ~ group, data= side_data)
anovamod <- oneway.test(drinksperday.diff  ~ group , data= side_data, var.equal = TRUE)

Fstat <- anovamod$statistic
J <- dim(Xmatrix)[2]
N <- dim(Xmatrix)[1]
Delta <- 0.01

pf(Fstat, J - 1, df2 = N-J, ncp = (Delta * N) / (1 - Delta), lower.tail=TRUE)
# or using our custom function:
noninf_ANOVA(side_data$drinksperday.diff, side_data$group, delta=0.01)

# Here is an analogous Bayesian testing scheme to consider:
complete_side_data <- na.omit(data.frame(side_data$drinksperday.diff, side_data$group)); colnames(complete_side_data)<-c("drinksperday.diff","group")

anovaBF(drinksperday.diff  ~ group , data= complete_side_data)

#############################################################
# Here is how the (recommended) the non-inferiority analysis
# in an ANOVA context with unequal variances (drinksperday):

Xmatrix <- model.matrix(drinksperday.diff  ~ group, data= side_data)
anovamod <- oneway.test(drinksperday.diff  ~ group , data= side_data, var.equal = FALSE)

Fprime <- anovamod$statistic
dfprime <-  anovamod$parameter[2]
J <- dim(Xmatrix)[2]
N <- dim(Xmatrix)[1]
Delta <- 0.01

pf(Fprime, J - 1, df2 = dfprime, ncp = (Delta * N) / (1 - Delta), lower.tail=TRUE)
# or using our custom function:
noninf_wellekANOVA(side_data$drinksperday.diff, side_data$group, delta=0.01)


#############################################################
###  Here is one alternative way of doing the non-inferiority analysis
# in a linear regression context (with 2 nested models) (drinksperday):

Xmatrix0 <- model.matrix(drinksperday.diff  ~ drinksperday.x, data= side_data)
Xmatrix1 <- model.matrix(drinksperday.diff  ~ drinksperday.x + group, data= side_data)

lmmodel0 <- lm(drinksperday.y  ~  drinksperday.x, data= side_data)
lmmodel1 <- lm(drinksperday.y  ~ drinksperday.x + group , data= side_data)

R2_diff <- summary(lmmodel1)$r.squared - summary(lmmodel0)$r.squared

Fstat <- anova(lmmodel1, lmmodel0)$F[2]
K <- dim(Xmatrix1)[2] - dim(Xmatrix0)[2] 
N <- dim(Xmatrix1)[1]
Delta <- 0.01
 
pf(Fstat,df1=K,df2=N-K-1,ncp=(N*Delta)/(1-Delta),lower.tail=TRUE)



