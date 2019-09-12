# source('~/Box Sync/Project/Seojin/mydistofy.R', chdir = TRUE)
###########################
## Set Working Directory ##
###########################

setwd("C:/Users/Seojin/Downloads")
setwd("~/Box Sync/Project/")
source('Seojin/myestimates.R')

##############
## Packages ##
##############
library(survival)
library(ggplot2)
library(OIsurv)
library(muhaz)


##########
## Data ##
##########

datp <- read.table("Data/wpbc.newdata.txt", header = T)	# 198 * 35
myidx <- read.table("Data/TenFoldCV_SampleIdx.txt", header = T)
datp1 <- datp[which(datp$outcome == 1),]
#datp <- read.table("wpbc.newdata.txt", header = T)	# 198 * 35
#myidx <- read.table("TenFoldCV_SampleIdx.txt", header = T)
#datp1 <- datp[which(datp$outcome == 1),]


#######################
## Distribution of Y ##
#######################

## Set Survival Data
my.surv <- Surv(datp$time, datp$outcome)

## Kaplan-Meier estimate for survival function (i.e. No-features)
km.fit <- survfit(my.surv~1)

## KM estimates for all time at "x"
mytime = 0:300
mysurvival <- myestimates(x = mytime, time = km.fit$time, estimates = km.fit$surv)
write.table(data.frame(time = mytime, survivalProb = mysurvival), file.path("Result/Distributions/SurvivalProbability.txt"), col.names = T, row.names = F)

## Hazard Function estimates for all time at "x"
mytime = 0:300
fit1 <- muhaz(datp$time, datp$outcome)
myhazard <- myestimates(x = mytime, time = fit1$est.grid, estimates = fit1$haz.est, type = "hazard")
write.table(data.frame(time = mytime, hazard = myhazard), file.path("Result/Distributions/HazardFunction.txt"), col.names = T, row.names = F)
ggplot(data.frame(time = fit1$est.grid, haz = fit1$haz.est), aes(x = time, y = haz)) + 
  geom_line() + ylab("Hazard Function")  


## plot: Distribution of Y with DFS
par(mfrow = c(5,5))
y = 1:130
for (dfs in seq(0, 120, by = 5)) {
	par(mar = c(1,1,3,1))
	a <- myestimates(x = dfs + y, time = fit1$est.grid, estimates = fit1$haz.est)
	b <- myestimates(x = dfs, time = km.fit$time, estimates = km.fit$surv)
	plot(y, a/b, xlab = "y", ylab = "f(y)", main = paste0("Dist of y under DFS ", dfs), type = 'l')
}

