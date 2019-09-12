# source("SurvivalAnalysis_seojin.R")
###########################
## Set Working Directory ##
###########################

setwd("~/Dropbox/0MachineLearning/Project_seojin/")
source("codes/ggsurv.R")
source("codes/myPHtest.R")
source("codes/mydrawcoxnell.R")
source("codes/mydrawmartingale.R")
source("codes/mydrawdeviance.R")
source("codes/mydrawschoenfeld.R")
source("codes/mydrawdfbeta.R")
source("codes/mylpauc.R")
source("codes/mymodelchecking.R")
source("codes/myerror.R")

##############
## Packages ##
##############

#install.packages("survival")
#install.packages("OIsurv")
#install.packages("MASS")
#install.packages("reshape")
#install.packages("survAUC")
#install.packages("muhaz")
library(survival)
library(OIsurv)
library(MASS)
library(reshape)
library(ggplot2)
library(muhaz)
library(survAUC)

##########
## Data ##
##########

datp <- read.table("Data/wpbc.newdata.txt", header = T)	# 198 * 35

###################
## Cox Model Fit ##
###################

## Set Survival Data
my.surv <- Surv(datp$time, datp$outcome)

## Kaplan-Meier estimate for survival function (i.e. No-features)
km.fit <- survfit(my.surv~1)
ggsurv(km.fit)

## Cox PH model
coxph.fit <- coxph(my.surv ~ . , method="breslow", data = datp[4:ncol(datp)])
coxph.fit.subAIC <- step(coxph.fit)	# AIC based feature selection
coxph.fit.subBIC <- step(coxph.fit, k = log(coxph.fit$n))	# BIC based feature selection
anova(coxph.fit.subAIC, coxph.fit) # chi-sq : 0.9619
anova(coxph.fit.subBIC, coxph.fit) # chi-sq : 0.1032
anova(coxph.fit.subBIC, coxph.fit.subAIC) # chi-sq : 0.001335

# CoxPH Full model
#Concordance= 0.791  (se = 0.045 )
#Rsquare= 0.264   (max possible= 0.899 )
#Likelihood ratio test= 60.81  on 32 df,   p=0.001573
#Wald test            = 38.91  on 32 df,   p=0.1867
#Score (logrank) test = 65.37  on 32 df,   p=0.0004501

# CoxPH sub model AIC
#Concordance= 0.768  (se = 0.045 )
#Rsquare= 0.228   (max possible= 0.899 )
#Likelihood ratio test= 51.19  on 13 df,   p=1.859e-06
#Wald test            = 50.09  on 13 df,   p=2.881e-06
#Score (logrank) test = 56.75  on 13 df,   p=1.991e-07

# CoxPH sub model BIC
#Concordance= 0.703  (se = 0.045 )
#Rsquare= 0.099   (max possible= 0.899 )
#Likelihood ratio test= 20.72  on 2 df,   p=3.167e-05
#Wald test            = 22.08  on 2 df,   p=1.602e-05
#Score (logrank) test = 23.45  on 2 df,   p=8.076e-06

## Check the Proportionality Assumption
my.PHtest(coxph.fit, data = datp) # perimeter3 is not ok
my.PHtest(coxph.fit.subAIC, data = datp) # every features is ok
my.PHtest(coxph.fit.subBIC, data = datp) # concavity3, GLOBAL is not ok

## Check Overall Model Fit via Residuals
model <- coxph.fit
model <- coxph.fit.subAIC
model <- coxph.fit.subBIC
#(a) Generalized (Cox-Snell) Residuals : check if x = y
my.draw.coxsnell(model, data = datp)
#(b) martingale -> check there is any other non-linear relation
my.draw.martingale(model, data = datp)
#(c) deviance
my.draw.deviance(model, data = datp)
#(d) Schoenfeld useful for assessing time trend or lack or proportionality, based on plotting versus event time
my.draw.schoenfeld(model)
#(e) dfbeta
my.draw.dfbeta(model)


################
## Prediction ##
################

## 10-Fold Training Data / Test Data Idx
#myidx <- data.frame(idx = sample(1:198, 198), 
#                  group = rep(1:10, each = 20)[1:nrow(datp)])
#myidx$ID <- datp$ID[myidx$idx]
#myidx <- myidx[order(myidx$idx),]
#write.table(myidx, "TenFoldCV_SampleIdx.txt", col.names = T, row.names = F)
myidx <- read.table("TenFoldCV_SampleIdx.txt", header = T)

## LP AUC
mylpauc.full <- my.lp.auc(data = datp, cvidx = myidx, 
                          modeltype = "Full", type = "cox")
apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
#sh         hc        Uno 
#0.7965011 11.5035097  0.6056163 
mylpauc.AIC <- my.lp.auc(data = datp, cvidx = myidx, modeltype = "AIC", type = "cox")
apply(mylpauc.AIC, 2, mean) # average auc over the 10-CV set
#sh         hc        Uno 
#0.7673026 12.0017438  0.6458688  
mylpauc.BIC <- my.lp.auc(data = datp, cvidx = myidx, modeltype = "BIC", type = "cox")
apply(mylpauc.BIC, 2, mean) # average auc over the 10-CV set
#sh        hc       Uno 
#0.7107065 9.7166019 0.5801592 

predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times,
        type = "brier", int.type = "unweighted")
predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times,
        type = "robust", int.type = "unweighted")

## Prediction
# library(muhaz)
# TR <- data[c(cvidx$idx[cvidx$group != i]),]
# TE <- data[c(cvidx$idx[cvidx$group == i]),]
# train.surv <- Surv(TR$time, TR$outcome)
# train.fit.all <- coxph(train.surv ~ . , method="breslow", 
#                        data = TR[,4:ncol(data)])
# train.fit <- step(train.fit.all)
#   
# my.pred.ch <- predict(train.fit, type = "expected")   # cumulative hazard: the expected number of events given the covariates and follow-up time
# my.pred.sp <- exp(-predict(train.fit, type = "expected"))   # estimated survival probability
# my.pred.lp <- predict(train.fit, newdata = TE)
# my.base.hazard <- muhaz(TR$time, TR$outcome, n.est.grid=112)
# my.hazard <- my.base.hazard$haz.est[match(TR$time, my.base.hazard$est.grid)]
# my.hazard[is.na(my.hazard)] <- 0
# my.pred.hazard <- sapply(exp(my.pred.lp), 
#                         f <- function(x, y) x*y, 
#                         y = my.hazard)
# my.pred.f <- apply(my.pred.hazard, 2, f <- function(x, y) x*y, y = my.pred.sp)
# my.pred.f <- apply(my.pred.f, 2, f <- function(x) x/sum(x))
# my.pred.time <- apply(my.pred.f, 2, f <- function(x, y) sum(x*y), y = TR$time)

##########################
## Parametric Model Fit ##
##########################
mytype = c("exponential", "weibull", "gaussian", "lognormal", "logistic", "loglogistic")

## Parametric model "weibull", "exponential", "gaussian", "logistic","lognormal" and "loglogistic".
fit0 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "exponential") #
fit1 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "weibull") #
fit2 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "gaussian")
fit3 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "lognormal") #
fit4 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "logistic")
fit5 <- survreg(my.surv ~ ., data = datp[,4:ncol(datp)], dist = "loglogistic") #
fit0.subAIC <- step(fit0)
fit1.subAIC <- step(fit1)
fit2.subAIC <- step(fit2)
fit3.subAIC <- step(fit3)
fit4.subAIC <- step(fit4)
fit5.subAIC <- step(fit5)
fit0.subBIC <- step(fit0, k = log(nrow(fit0$y)))
fit1.subBIC <- step(fit1, k = log(nrow(fit1$y)))
fit2.subBIC <- step(fit2, k = log(nrow(fit2$y)))
fit3.subBIC <- step(fit3, k = log(nrow(fit3$y)))
fit4.subBIC <- step(fit4, k = log(nrow(fit4$y)))
fit5.subBIC <- step(fit5, k = log(nrow(fit5$y)))

## get AIC, BIC
for (num in c(0:5)) {
	eval(parse(text = paste("fit <- fit", num, sep = "")))
	eval(parse(text = paste("fit.subAIC <- fit", num, ".subAIC", sep = "")))
	eval(parse(text = paste("fit.subBIC <- fit", num, ".subBIC", sep = "")))
	res.aic <- c(extractAIC(fit)[2], 
	             extractAIC(fit.subAIC)[2], 
	             extractAIC(fit.subBIC)[2])
	res.bic <- c(extractAIC(fit, k = log(198))[2], 
	          extractAIC(fit.subAIC, k = log(198))[2], 
	          extractAIC(fit.subBIC, k = log(198))[2])
	resres <- rbind(res.aic, res.bic)
	colnames(resres) <- c("all", "subAIC", "subBIC")
	rownames(resres) <- c("AICvalues", "BICvalues")
	print(resres)
# Don't have any preference btw AIC and BIC. The number of features are not that different
}
               # all   subAIC   subBIC
# AICvalues 611.7686 558.2640 571.0517
# BICvalues 720.2814 601.0115 584.2048
               # all   subAIC   subBIC
# AICvalues 613.3088 560.2136 572.5202
# BICvalues 725.1098 606.2493 588.9615
               # all   subAIC   subBIC
# AICvalues 674.4273 625.6853 625.3709
# BICvalues 786.2284 684.8741 658.2536
               # all   subAIC   subBIC
# AICvalues 608.5610 562.9376 566.1995
# BICvalues 720.3621 618.8381 585.9291
               # all   subAIC   subBIC
# AICvalues 680.2220 628.7767 628.0676
# BICvalues 792.0231 687.9656 660.9502
               # all   subAIC   subBIC
# AICvalues 611.1858 561.8094 562.8367
# BICvalues 722.9869 617.7100 595.7193

## Goodnees of Fit test
anova(fit0.subBIC, fit0.subAIC)$Pr[2] # 0.0003219087
anova(fit1.subBIC, fit1.subAIC)$Pr[2] # 0.05373154
anova(fit2.subBIC, fit2.subAIC)$Pr[2] # 0.0003889984
anova(fit3.subBIC, fit3.subAIC)$Pr[2] # 0.04710742
anova(fit4.subBIC, fit4.subAIC)$Pr[2] # 0.05373154
anova(fit5.subBIC, fit5.subAIC)$Pr[2] # 0.035652
# all parametric model, AIC is better than BIC model.

## Parametic Assumption Check : Weibull, Log-logistic, Log-Normal is good, exponential is bad
model=survfit(coxph.fit,type='aalen')
mymodelchecking(model)
vec <- residuals(fit2)
y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
ggplot(data.frame(y = vec), aes(sample = y)) + 
  stat_qq() + geom_abline(slope = slope, intercept = int, colour = "red") +
  ggtitle("Gaussian Model Check")
#  Exponential   Weibull Log-Logistic Log-Normal
#r2   0.9460781 0.9741865    0.9772381  0.9809825

################
## Prediction ##
################

#cox, exponential, weibull, gaussian, lognormal or logistic?
myidx <- read.table("Data/TenFoldCV_SampleIdx.txt", header = T)

## LP AUC
myres <- NULL
for (ty in c("exponential", "weibull", "gaussian", "lognormal", "logistic", "loglogistic")) {
  mylpauc.full <- my.lp.auc(data = datp, cvidx = myidx, modeltype = "Full", type = ty)
  mylpauc.AIC <- my.lp.auc(data = datp, cvidx = myidx, modeltype = "AIC", type = ty)
  mylpauc.BIC <- my.lp.auc(data = datp, cvidx = myidx, modeltype = "BIC", type = ty)
  r1 <- apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
  r2 <- apply(mylpauc.AIC, 2, mean) # average auc over the 10-CV set
  r3 <- apply(mylpauc.BIC, 2, mean) # average auc over the 10-CV set
  res <- data.frame(Full = r1, AIC = r2, BIC = r3, type = ty)
  myres <- rbind(myres, res)
}
 
## Prediction Error
myerror <- NULL
mymax = Inf
mymax = 130
for (ty in c("exponential", "weibull", "lognormal", "loglogistic")) {
  mylpauc.full <- my.error(data = datp, cvidx = myidx, modeltype = "Full", type = ty, mymax = mymax)
  mylpauc.AIC <- my.error(data = datp, cvidx = myidx, modeltype = "AIC", type = ty, mymax = mymax)
  mylpauc.BIC <- my.error(data = datp, cvidx = myidx, modeltype = "BIC", type = ty, mymax = mymax)
  r1 <- apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
  r2 <- apply(mylpauc.AIC, 2, mean) # average auc over the 10-CV set
  r3 <- apply(mylpauc.BIC, 2, mean) # average auc over the 10-CV set
  res <- data.frame(Full = r1, AIC = r2, BIC = r3, type = ty)
  myerror <- rbind(myerror, res)
}

#> myres
# Full       AIC BIC        type
# sh          NaN       NaN 0.5 exponential
# hc   10.0910524 9.8801931 0.0 exponential
# Uno   0.4049274 0.3762688 0.5 exponential
# sh1         NaN       NaN 0.5     weibull
# hc1  10.1673829 9.9126030 0.0     weibull
# Uno1  0.4150109 0.3733772 0.5     weibull
# sh2   0.9725011 0.9746276 0.5    gaussian
# hc2  10.6058150 9.8381625 0.0    gaussian
# Uno2  0.3726739 0.3354190 0.5    gaussian
# sh3         NaN       NaN 0.5   lognormal
# hc3  10.5691545 9.7984965 0.0   lognormal
# Uno3  0.3805698 0.3977294 0.5   lognormal
# sh4   0.9705041 0.9745181 0.5    logistic
# hc4  11.0842711 8.2857907 0.0    logistic
# Uno4  0.3875824 0.3635284 0.5    logistic
# sh        NaN        NaN 0.5 loglogistic
# hc  9.7277786 10.9657271 0.0 loglogistic
# Uno 0.3974762  0.4082817 0.5 loglogistic

#> myerror (old version error)
                # Full         AIC       BIC        type
# All       19.8748329  28.1953513 10.215149 exponential
# NonRecur   1.5461096   1.7865055  0.000000 exponential
# Recur     76.2869672 104.1332240 41.140184 exponential
# All1      22.6825956  24.8075975 14.979819     weibull
# NonRecur1  1.5069299   1.7917093  0.000000     weibull
# Recur1    88.2035827  92.1371632 60.270748     weibull
# All2      14.3670078  21.7358426 11.919007   lognormal
# NonRecur2  2.3831016   2.5543975  0.000000   lognormal
# Recur2    52.2581851  77.4444515 47.817603   lognormal
# All3      14.0971911  20.9397058  9.982416 loglogistic
# NonRecur3  2.6855415   2.6736252  0.000000 loglogistic
# Recur3    49.7360999  73.6075008 40.106666 loglogistic

#> myerror (new version error): some mistakes in the code
#               Full        AIC       BIC        type
#All       19.884039  28.210873 10.215149 exponential
#NonRecur   1.546110   1.786505  0.000000 exponential
#Recur     76.332997 104.210833 41.140184 exponential
#All1      22.699122  24.823508 14.979819     weibull
#NonRecur1  1.506930   1.791709  0.000000     weibull
#Recur1    88.286212  92.216717 60.270748     weibull
#All2      14.412804  21.274401 11.919007   lognormal
#NonRecur2  2.383102   2.621889  0.000000   lognormal
#Recur2    52.479026  75.623182 47.817603   lognormal
#All3      14.152535  20.979935  9.982416 loglogistic
#NonRecur3  2.685541   2.673625  0.000000 loglogistic
#Recur3    50.002575  73.791655 40.106666 loglogistic

#> myerror (code error fixed version. no maximum predicted time)
#                Full        AIC       BIC        type
#All        77.879107 105.997338  41.14018 exponential
#NonRecur    2.022392   2.282406   0.00000 exponential
#Recur     325.207150 419.542895 171.91732 exponential
#All1       89.793142  94.008426  60.27075     weibull
#NonRecur1   1.974162   2.289250   0.00000     weibull
#Recur1    377.197823 376.989912 251.54831     weibull
#All2       54.862128  80.096248  47.81760   lognormal
#NonRecur2   3.144198   3.279324   0.00000   lognormal
#Recur2    231.498856 314.433687 198.99551   lognormal
#All3       52.688117  76.465280  40.10667 loglogistic
#NonRecur3   3.507189   3.452024   0.00000 loglogistic
#Recur3    217.557650 297.311326 167.14682 loglogistic

#> myerror (code error fixed version. by forcing maximum predicted time to 130)
#                Full        AIC       BIC        type
#All       21.700601 20.215595  24.96889 exponential
#NonRecur   2.022392  2.282406   0.00000 exponential
#Recur     84.185525 77.462583 104.93833 exponential
#All1      21.975712 20.832120  24.96889     weibull
#NonRecur1  1.974162  2.289250   0.00000     weibull
#Recur1    85.461880 79.703218 104.93833     weibull
#All2      21.037779 20.678357  24.96889   lognormal
#NonRecur2  3.144198  3.279324   0.00000   lognormal
#Recur2    78.044954 75.283070 104.93833   lognormal
#All3      20.773594 19.557755  24.96889 loglogistic
#NonRecur3  3.507189  3.452024   0.00000 loglogistic
#Recur3    75.768623 70.882762 104.93833 loglogistic
####################################################

## fit0
## AIC/BIC of fit, fit.subAIC, fit.subBIC
#[1] 611.7686 558.2640 571.0517
#[1] 720.2814 601.0115 584.2048
my.surv ~ radius1 + radius3 + texture2 + perimeter1 + 
  area1 + area2 + area3 + compactness1 + compactness2 + symmetry2 + 
  symmetry3 + lymphNodeStatus
my.surv ~ radius1 + radius3 + area1

## fit1
# AIC of fit, fit.subAIC, fit.subBIC
# 613.3088 560.2136 572.5202
# BIC of fit, fit.subAIC, fit.subBIC
# 725.1098 606.2493 588.9615
Step:  AIC=560.21
my.surv ~ radius1 + radius3 + texture2 + perimeter1 + area1 + 
    area2 + area3 + compactness1 + compactness2 + symmetry2 + 
    symmetry3 + lymphNodeStatus
Step:  BIC=588.9615
my.surv ~ radius1 + radius3 + area1

## fit2
## AIC/BIC of fit, fit.subAIC, fit.subBIC
#[1] 674.4273 625.6853 625.3709
#[1] 786.2284 684.8741 658.2536
Step:  AIC=625.69
my.surv ~ radius1 + radius3 + perimeter1 + area1 + area2 + area3 + 
    smoothness2 + compactness1 + compactness2 + compactness3 + 
    concavity1 + symmetry1 + symmetry2 + symmetry3 + fractalDimension2 + 
    lymphNodeStatus
Step:  BIC=658.2536
my.surv ~ area1 + area2 + compactness1 + compactness2 + 
  symmetry1 + symmetry2 + symmetry3 + lymphNodeStatus

## fit3
## AIC/BIC of fit, fit.subAIC, fit.subBIC
#[1] 608.5610 562.9376 566.1995
#[1] 720.3621 618.8381 585.9291
Step:  AIC=562.94
my.surv ~ radius1 + radius3 + perimeter3 + area1 + area2 + area3 + 
    smoothness2 + compactness1 + compactness2 + compactness3 + 
    concavity1 + symmetry1 + symmetry2 + symmetry3 + lymphNodeStatus
Step:  AIC=585.9291
my.surv ~ area1 + area2 + symmetry1 + lymphNodeStatus

## fit4
## AIC/BIC of fit, fit.subAIC, fit.subBIC
#[1] 680.2220 628.7767 628.0676
#[1] 792.0231 687.9656 660.9502
Step:  AIC=628.78
my.surv ~ radius1 + radius3 + perimeter1 + area1 + area2 + area3 + 
    smoothness2 + compactness1 + compactness2 + compactness3 + 
    concavity1 + symmetry1 + symmetry2 + symmetry3 + fractalDimension2 + 
    lymphNodeStatus
Step:  AIC=660.9502
my.surv ~ area1 + area2 + compactness1 + compactness2 + 
  symmetry1 + symmetry2 + symmetry3 + lymphNodeStatus

## fit5
## AIC/BIC of fit, fit.subAIC, fit.subBIC
# AICvalues 611.1858 561.8094 562.8367
# BICvalues 722.9869 617.7100 595.7193
Step:  AIC=561.8094
my.surv ~ radius1 + radius3 + perimeter3 + 
    area1 + area2 + area3 + smoothness2 + compactness1 + compactness2 + 
    compactness3 + concavity1 + symmetry1 + symmetry2 + symmetry3 + 
    lymphNodeStatus








