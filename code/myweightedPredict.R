#source("myweightedPredict.R")

myweightedPredict0 <- function(beta, test, type = "gaussian", modeltype = "Full", mymax = 130) {
	mybeta <- as.numeric(beta)
	mytest <- cbind(1, test)
	xb <- as.matrix(mytest) %*% as.matrix(mybeta)

	if (type == "gaussian") {
		tnew <- exp(xb)
	} else if (type == "gamma_log") {
		tnew <- exp(xb)
	} else {
		stop("Put your model type : gaussian, or gamma_log ?")
	}
	tnew <- sapply(tnew, min, mymax)
}

myweightedPredict <- function(beta, test, type = "gaussian", mymax = 130, modeltype = "NoVar") {
	myvar <- names(beta)[-1]
	mybeta <- as.numeric(beta)
	if (modeltype == "Full") {
		mytest <- cbind(1, test)
		xb <- as.matrix(mytest) %*% as.matrix(mybeta)
	} else if (modeltype == "AIC") {
		myvar <- gsub("as.matrix\\(TR\\)", "", myvar)
		mytest <- test[,myvar]
		mytest <- cbind(1, mytest)
		xb <- as.matrix(mytest) %*% as.matrix(mybeta)
	} else if (modeltype == "BIC") {
		myvar <- gsub("as.matrix\\(TR\\)", "", myvar)
		mytest <- test[,myvar]
		mytest <- cbind(1, mytest)
		xb <- as.matrix(mytest) %*% as.matrix(mybeta)
	} else {
		stop("Put your model type : Full, AIC, or BIC")
	}

	if (type == "gaussian") {
		tnew <- exp(xb)
	} else if (type == "gamma_log") {
		tnew <- exp(xb)
	} else {
		stop("Put your model type : gaussian, or gamma_log ?")
	}
	#tnew <- sapply(tnew, min, mymax)
}