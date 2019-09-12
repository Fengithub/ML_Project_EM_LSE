# source("OriginalLinearModel_seojin.R")
###########################
## Set Working Directory ##
###########################

setwd("~/Box Sync/Project/")
source("Seojin/mylinearerror.R")

##############
## Packages ##
##############


##########
## Data ##
##########

datp <- read.table("Data/wpbc.newdata.txt", header = T)	# 198 * 35
myidx <- read.table("Data/TenFoldCV_SampleIdx.txt", header = T)
datp1 <- datp[which(datp$outcome == 1),]

##BIC 
varselected <- c("radius2", "radius3", "texture2", "perimeter1", "perimeter2", "perimeter3", "area1", "area2", "smoothness1", "smoothness2", "compactness2", "compactness3",
					"concavity1", "concavePoints1", "concavePoints2", "concavePoints3", "symmetry1", "fractalDimension2", "fractalDimension3")
datp <- datp[,c("ID", "outcome", "time", varselected)]


######################
## Linear Model Fit ##
######################

lin.fit0 <- lm(log(time) ~ 1, data = datp1[,-c(1,2)])
lin.fit <- lm(log(time) ~ .,data = datp1[,-c(1,2)])
lin.fit.subAIC <- step(lin.fit)
lin.fit.subBIC <- step(lin.fit, k = log(length(lin.fit$residuals)))
summary(lin.fit0)
summary(lin.fit)
summary(lin.fit.subAIC)
summary(lin.fit.subBIC)
# AIC model:
#log(time) ~ radius2 + radius3 + texture2 + texture3 + 
#    perimeter1 + perimeter2 + perimeter3 + area1 + area2 + area3 + 
#    smoothness1 + smoothness2 + compactness2 + compactness3 + 
#    concavity1 + concavePoints1 + concavePoints2 + concavePoints3 + 
#    symmetry1 + fractalDimension2 + fractalDimension3
#BIC model:
#log(time) ~ radius2 + radius3 + texture2 + perimeter1 + 
#    perimeter2 + perimeter3 + area1 + area2 + smoothness1 + smoothness2 + 
#    compactness2 + compactness3 + concavity1 + concavePoints1 + 
#    concavePoints2 + concavePoints3 + symmetry1 + fractalDimension2 + 
#    fractalDimension3

lin.loggamma.fit0 <- glm(time ~ 1,data = datp1[,-c(1,2)], family = Gamma(link = "log")) 
lin.loggamma.fit <- glm(time ~ .,data = datp1[,-c(1,2)], family = Gamma(link = "log"))
lin.loggamma.fit.subAIC <- step(lin.loggamma.fit)
lin.loggamma.fit.subBIC <- step(lin.loggamma.fit, k = log(length(lin.loggamma.fit$residuals)))
summary(lin.loggamma.fit)
summary(lin.loggamma.fit.subAIC)
summary(lin.loggamma.fit.subBIC)
#AIC model:
#radius2 + radius3 + texture2 + texture3 + 
#    perimeter1 + perimeter2 + perimeter3 + area1 + area2 + smoothness1 + 
#    smoothness2 + compactness2 + compactness3 + concavity1 + 
#    concavePoints1 + concavePoints2 + concavePoints3 + symmetry1 + 
#    fractalDimension2 + fractalDimension3
#BIC model:
#radius2 + radius3 + texture2 + perimeter1 + 
#    perimeter2 + perimeter3 + area1 + area2 + smoothness1 + smoothness2 + 
#    compactness2 + compactness3 + concavity1 + concavePoints1 + 
#    concavePoints2 + concavePoints3 + symmetry1 + fractalDimension2 + 
#    fractalDimension3

## Model Assumption Check
par(mfrow = c(2,2))
plot(lin.fit)
plot(lin.fit.subAIC)
plot(lin.fit.subBIC)
plot(lin.loggamma.fit)
plot(lin.loggamma.fit.subAIC)
plot(lin.loggamma.fit.subBIC)

## Goodness of Fit Test
anova(lin.fit.subAIC, lin.fit) # 0.9952 : AIC preferred
anova(lin.fit.subBIC, lin.fit) # 0.9827 : BIC preferred
anova(lin.fit.subBIC, lin.fit.subAIC) # 0.2833: BIC preferred
a1 <- anova(lin.loggamma.fit.subAIC, lin.loggamma.fit); 1-pchisq(a1$Deviance[2], df = a1$Df[2])# 0.9997183 : AIC preferred
a1 <- anova(lin.loggamma.fit.subBIC, lin.loggamma.fit); 1-pchisq(a1$Deviance[2], df = a1$Df[2]) # 0.9996521 : BIC preferred
a1 <- anova(lin.loggamma.fit.subBIC, lin.loggamma.fit.subAIC); 1-pchisq(a1$Deviance[2], df = a1$Df[2]) # 0.5140017: BIC preferred

################
## Prediction ##
################
source("Seojin/mylinearerror.R")

mymax = 130
myerror <- NULL
for (ty in c("gaussian", "gamma_log")) {
  mylpauc.null <- my.linearerror(data = datp, cvidx = myidx, modeltype = "NoVar", type = ty, mymax = mymax)
  mylpauc.full <- my.linearerror(data = datp, cvidx = myidx, modeltype = "Full", type = ty, mymax = mymax)
  mylpauc.AIC <- my.linearerror(data = datp, cvidx = myidx, modeltype = "AIC", type = ty, mymax = mymax)
  mylpauc.BIC <- my.linearerror(data = datp, cvidx = myidx, modeltype = "BIC", type = ty, mymax = mymax)
  r0 <- apply(mylpauc.null, 2, mean)
  r1 <- apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
  r2 <- apply(mylpauc.AIC, 2, mean) # average auc over the 10-CV set
  r3 <- apply(mylpauc.BIC, 2, mean) # average auc over the 10-CV set
  res <- data.frame(Novar = r0, Full = r1, AIC = r2, BIC = r3, type = ty)
  myerror <- rbind(myerror, res)
}

mymax = 130
myerror <- NULL
for (ty in c("gaussian", "gamma_log")) {
  mylpauc.null <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "NoVar", type = ty, mymax = mymax)
  mylpauc.full <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "Full", type = ty, mymax = mymax)
  r0 <- apply(mylpauc.null, 2, mean)
  r1 <- apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
  res <- data.frame(Novar = r0, Full = r1, type = ty)
  myerror <- rbind(myerror, res)
}

mymax = 130
myerror <- NULL
for (ty in c("gaussian", "gamma_log")) {
  mylpauc.null <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "NoVar", type = ty, mymax = mymax)
  mylpauc.full <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "Full", type = ty, mymax = mymax)
  mylpauc.AIC <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "AIC", type = ty, mymax = mymax)
  mylpauc.BIC <- my.linearerror.final(data = datp, cvidx = myidx, modeltype = "BIC", type = ty, mymax = mymax)
  r0 <- apply(mylpauc.null, 2, mean)
  r1 <- apply(mylpauc.full, 2, mean) # average auc over the 10-CV set
  r2 <- apply(mylpauc.AIC, 2, mean) # average auc over the 10-CV set
  r3 <- apply(mylpauc.BIC, 2, mean) # average auc over the 10-CV set
  res <- data.frame(Novar = r0, Full = r1, AIC = r2, BIC = r3, type = ty)
  myerror <- rbind(myerror, res)
}

for (ty in c("gaussian", "gamma_log")) {
	#mypred.null <- my.linearpredicted(data = datp, cvidx = myidx, modeltype = "NoVar", type = ty)
	mypred.full <- my.linearpredicted(data = datp, cvidx = myidx, modeltype = "Full", type = ty)
	#mypred.AIC <- my.linearpredicted(data = datp, cvidx = myidx, modeltype = "AIC", type = ty)
	#mypred.BIC <- my.linearpredicted(data = datp, cvidx = myidx, modeltype = "BIC", type = ty)
}

write.table(mypred.full, "Result/Baseline_Linear(Original Model)/GaussianLinearModel_iter1_BIC.txt")

#> myerror (This has some mistakes)
#               Full       AIC       BIC      type
#All       31.122489 31.093246 30.893729  gaussian
#NonRecur  38.956518 38.918149 38.620542  gaussian
#Recur      5.269092  5.271378  5.330927  gaussian
#All1      31.133423 31.091039 30.602678 gamma_log
#NonRecur1 38.974734 38.922792 38.271185 gamma_log
#Recur1     5.256083  5.241781  5.234937 gamma_log

#> myerror (no maximum predicted time)
#                 Full        AIC      BIC      type
#All        1480.41523  2015.8617 33.80984  gaussian
#NonRecur     33.67772    33.2156 38.94909  gaussian
#Recur      5839.06466  7986.6153 17.02149  gaussian
#All1       3208.65127  4977.7254 28.77371 gamma_log
#NonRecur1    33.56749    33.3144 31.81935 gamma_log
#Recur1    12757.77246 19837.4243 18.52445 gamma_log

#> myerror (By forcing maximum predicted time to be 130)
#              Full      AIC      BIC      type
#All       33.66711 32.77410 33.80984  gaussian
#NonRecur  33.67772 33.21560 38.94909  gaussian
#Recur     35.92717 33.15641 17.02149  gaussian
#All1      33.72715 33.49348 28.77371 gamma_log
#NonRecur1 33.56749 33.31440 31.81935 gamma_log
#Recur1    36.49882 36.43061 18.52445 gamma_log

##
source("Seojin/mylinearerror.R")
mypred <- NULL
for (i in 1:10) {
	res <- my.linearerror.sub.final(i, data = datp, cvidx = myidx, modeltype = "Full", type = ty, mymax = 130)
	mypred <- rbind(mypred, res)
}
write.table(mypred, "Result/Baseline_Linear(Original Model)/GaussianLinearModel_iter1_nomax.txt", row.names = F, col.names = T, quote = F)

## error
source("Seojin/myerror_general.R")
res <- read.table("Result/Baseline_Linear(Original Model)/GaussianLinearModel_iter1.txt", header = T)
res <- read.table("Result/Baseline_Linear(Original Model)/GaussianLinearModel_iter1_nomax.txt", header = T)
myerror <- myerrorfinal(predicted = res$predicted, outcome = res$outcome, time = res$time)
myerror
exp(myerror) # all nonrecur recur 

###### error
## max 130
> myerror
[1] 0.2026832 0.2295492 0.1163688
> exp(sqrt(myerror)) # all nonrecur recur
[1] 1.568631 1.614642 1.406534

#no max
> myerror
[1] 0.2044819 0.2310106 0.1192515
> exp(sqrt(myerror)) # all nonrecur recur 
[1] 1.571761 1.617102 1.412453


###### check mannually for no max
dfs.observed = res$time[outcome == 0] # original scale
ttr.predicted = res$predicted[outcome == 0] # original scale

dfs.observed = 9
ttr.predicted = 2.702613e+02
mysurvival.list <- read.table("Result/Distributions/SurvivalProbability.txt", header = T)
myhazard.list <- read.table("Result/Distributions/HazardFunction.txt", header = T)
y = 1:130








