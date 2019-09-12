#source("EM_iterative.R")
###########################
## Set Working Directory ##
###########################

setwd("~/Box Sync/Project/")
#setwd("C:/Users/Seojin/Downloads/")
#source("Seojin/myweightedLSE.R")
source("Seojin/myweightedPredict.R")
source("Seojin/myerror_general.R")
source("Seojin/myEMiter.R")

##########
## Data ##
##########

datp <- read.table("Data/wpbc.newdata.txt", header = T)	# 198 * 35
cvidx <- read.table("Data/TenFoldCV_SampleIdx.txt", header = T)
#datp <- read.table("wpbc.newdata.txt", header = T)	# 198 * 35
#cvidx <- read.table("TenFoldCV_SampleIdx.txt", header = T)

##BIC 
varselected <- c("radius2", "radius3", "texture2", "perimeter1", "perimeter2", "perimeter3", "area1", "area2", "smoothness1", "smoothness2", "compactness2", "compactness3",
					"concavity1", "concavePoints1", "concavePoints2", "concavePoints3", "symmetry1", "fractalDimension2", "fractalDimension3")
datp <- datp[,c("ID", "outcome", "time", varselected)]


########
## EM ##
########
dat = datp
thres = 10^-2
niter = 5
mymodeltype = "Full"
mymodeltype0 = "BIC"

	eval(parse(text = paste0("mypred.final.", mymodeltype," <- NULL")))
	eval(parse(text = paste0("mypred.final.train", mymodeltype," <- NULL")))
	for (i in 1:10) {
		cat(i, "\n")
  		TR <- datp[c(cvidx$idx[cvidx$group != i]),]
  		TE <- datp[c(cvidx$idx[cvidx$group == i]),]
		eval(parse(text = paste0("res1 <- read.table('Result/WeightedLSE/WeightedLSE_iter1_", mymodeltype0,"_trainCV",i,".txt', header = T)")))
		nextbeta = 0
		myres <- NULL
		myres$beta <- 100
		for(j in 1:niter){
			myres <- myEMiter(time = res1$observedtime, outcome = res1$observedoutcome, predicted = res1$predictedtime, 
						id = res1$indexforpatiemt , modeltype = mymodeltype, type = "gaussian")
			mypred <- myweightedPredict(beta = myres$beta, test = TE[,4:ncol(TE)], type = "gaussian", mymax = 130, 
						modeltype = mymodeltype)
			mypred.train <- myweightedPredict(beta = myres$beta,
						test = TR[,4:ncol(TR)],
						type = "gaussian", modeltype = mymodeltype) 
			myerror <- myerrorfinal(predicted = mypred, outcome = TE$outcome, time = TE$time)
			cat("error:", myerror, "\n")
			cat("scaled error:", exp(sqrt(myerror)), "\n")
			if ( j == niter) {
				#write.table("")
				myid.final <- TE[,c("ID", "outcome", "time")]
				myid.final$predicted <- mypred
				myid.final$group <- i
				myid.final.train <- TR[,c("ID", "outcome", "time")]
				myid.final.train$predicted <- mypred.train
				myid.final.train$group <- i
				eval(parse(text = paste0("mypred.final.", mymodeltype," <- rbind(mypred.final.", mymodeltype,", myid.final)")))
				eval(parse(text = paste0("write.table(myid.final.train, 'Result/WeightedLSE/WeightedLSE_", mymodeltype0,"trainCV", i,".txt', col.names = T, row.names = F, quote = F)")))
				eval(parse(text = paste0("write.table(myres$beta, 'Result/WeightedLSE/WeightedLSE_", mymodeltype0,"_betaCV", i,".txt', col.names = T, row.names = T, quote = F)")))
				#return(myid.final)
			}
			nextbeta = myres$beta
		}
	}

		final.res <- mypred.final.Full
		final.res <- final.res[,c('ID', 'group', 'outcome', 'time', 'predicted')]
		colnames(final.res) <- c('indexforpatiemt', 'idxgroup', 'observedoutcome', 'observedtime', 'predictedtime')
		eval(parse(text = paste0("write.table(final.res, 'Result/WeightedLSE/WeightedLSE_",mymodeltype0,".txt', col.names = T, row.names = F, quote = F)")))


a <- mypred.final.Full
myerror <- myerrorfinal(predicted = a$predicted, outcome = a$outcome, time = a$time)
myerror # all nonrecur recur 
exp(sqrt(myerror))

