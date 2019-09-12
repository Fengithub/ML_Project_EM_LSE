#source("EM_WeightLSE.R")
###########################
## Set Working Directory ##
###########################

setwd("~/Box Sync/Project/")
#setwd("C:/Users/Seojin/Downloads/")
source("Seojin/myweightedLSE.R")
source("Seojin/myweightedPredict.R")
source("Seojin/myerror_general.R")

##########
## Data ##
##########

datp <- read.table("Data/wpbc.newdata.txt", header = T)	# 198 * 35
cvidx <- read.table("Data/TenFoldCV_SampleIdx.txt", header = T)
haz <- read.table("Result/Distributions/HazardFunction.txt", header = T)
surv <- read.table("Result/Distributions/SurvivalProbability.txt", header = T)
#datp <- read.table("wpbc.newdata.txt", header = T)	# 198 * 35
#cvidx <- read.table("TenFoldCV_SampleIdx.txt", header = T)

##BIC 
varselected <- c("radius2", "radius3", "texture2", "perimeter1", "perimeter2", "perimeter3", "area1", "area2", "smoothness1", "smoothness2", "compactness2", "compactness3",
					"concavity1", "concavePoints1", "concavePoints2", "concavePoints3", "symmetry1", "fractalDimension2", "fractalDimension3")
datp <- datp[,c("ID", "outcome", "time", varselected)]

##################
## Weighted LSE ##
##################
mymodeltype = "Full"
mymodeltype0 = mymodeltype
mymodeltype0 = "BIC"
for (mymodeltype in c("AIC", "BIC")) {
	eval(parse(text = paste0("mypred.final.", mymodeltype," <- NULL")))
	eval(parse(text = paste0("mypred.final.train", mymodeltype," <- NULL")))
	for (i in 1:10) {
		cat(i, "\n")
  		TR <- datp[c(cvidx$idx[cvidx$group != i]),]
  		TE <- datp[c(cvidx$idx[cvidx$group == i]),]
		myfit <- myweightedLSE(time = TR$time, outcome = TR$outcome,
			train = TR[,4:ncol(TR)], 
			hazard = haz, survival = surv, y = 1:130,
			type = "gaussian", modeltype = mymodeltype)
		mybeta <- myfit$beta
		if (length(mybeta) > 1) {
			names(mybeta)[2:length(mybeta)] <- gsub("as.matrix\\(TR\\)", "", names(mybeta)[2:length(mybeta)])
			cat("myweightedLSE\n")
			mypred <- myweightedPredict(beta = myfit$beta,
				test = TE[,4:ncol(TE)],
				type = "gaussian", modeltype = mymodeltype)
			mypred.train <- myweightedPredict(beta = myfit$beta,
				test = TR[,4:ncol(TR)],
				type = "gaussian", modeltype = mymodeltype)
		} else {
			mypred = rep(mybeta, nrow(TE))
			mypred.train = rep(mybeta, nrow(TR))
		}
		myid.final <- TE[,c("ID", "outcome", "time")]
		myid.final$predicted <- mypred
		myid.final$group <- i
		myid.final.train <- TR[,c("ID", "outcome", "time")]
		myid.final.train$predicted <- mypred.train
		myid.final.train$group <- i
		myid.final.train <- myid.final.train[,c("ID", "group", "outcome", "time", "predicted")]
		colnames(myid.final.train) <- c("indexforpatiemt", "idxgroup", "observedoutcome", "observedtime", "predictedtime")
		eval(parse(text = paste0("mypred.final.", mymodeltype," <- rbind(mypred.final.", mymodeltype,", myid.final)")))
		eval(parse(text = paste0("write.table(myid.final.train, 'Result/WeightedLSE/WeightedLSE_iter1_",mymodeltype0,"_trainCV", i,".txt', col.names = T, row.names = F, quote = F)")))
		eval(parse(text = paste0("write.table(mybeta, 'Result/WeightedLSE/WeightedLSE_iter1_",mymodeltype0,"_betaCV", i,".txt', col.names = T, row.names = T, quote = F)")))
	}
	eval(parse(text = paste0("mypred.final.Full <- mypred.final.", mymodeltype,"[,c('ID', 'group', 'outcome', 'time', 'predicted')]")))
	eval(parse(text = paste0("colnames(mypred.final.", mymodeltype,") <- c('indexforpatiemt', 'idxgroup', 'observedoutcome', 'observedtime', 'predictedtime')")))
	eval(parse(text = paste0("write.table(mypred.final.", mymodeltype,", 'Result/WeightedLSE/WeightedLSE_iter1_",mymodeltype0,".txt', col.names = T, row.names = F, quote = F)")))
}

a <- mypred.final.BIC
a <- mypred.final.AIC
a <- mypred.final.Full
myerror <- myerrorfinal(predicted = a$predictedtime, outcome = a$observedoutcome, time = a$observedtime)
myerror # all nonrecur recur 
exp(sqrt(myerror))
#> myerror # all nonrecur recur 
#[1] 0.1633500 0.1578302 0.1810839
#> exp(sqrt(myerror))
#[1] 1.498052 1.487770 1.530416
#BIC
#> myerror # all nonrecur recur 
#[1] 0.1615224 0.1721903 0.1272488
#> exp(sqrt(myerror))
#[1] 1.494660 1.514308 1.428635


