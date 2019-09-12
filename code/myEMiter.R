#source("myEMiter.R")

myEMiter <- function(time, outcome, predicted, id, modeltype = "Full", type = "gaussian"){
	newy = predicted
	newy[outcome == 0] <- mapply(max, newy[outcome == 0], time[outcome == 0])
	mydat = dat[match(id , dat$ID),]

	if (type == "gaussian") {
		my.fit <- lm(log(newy) ~ ., data = mydat[,4:ncol(mydat)])
	} else if (type == "gamma_log") {
		my.fit <- glm(newy ~ ., data = mydat[,4:ncol(mydat)], family = Gamma(link = "log"))
	} else {
		stop("Put your model type : gaussian, or gamma_log?")
	}

	if (modeltype == "NoVar") {
    		train.fit <- NULL
		coef = "NULL"
  	} else if (modeltype == "AIC") {
    		train.fit <- step(my.fit)
		coefnew = names(train.fit$coefficients)
 	} else if (modeltype == "BIC") {
  		train.fit <- step(my.fit, k = nrow(TR))
		coefnew = names(train.fit$coefficients)
 	} else if (modeltype == "Full") {
  		train.fit <- my.fit
		coefnew = names(train.fit$coefficients)
 	} else {
		stop("Put your model type : AIC or BIC or Full or NoVar?")
 	}
 	myfitted = exp(train.fit$fitted.values)
 	myfitted = sapply(myfitted, min, 130)

	res = list(beta = train.fit$coefficients, coefficients = coefnew, fittedvalues = myfitted)
	return(res)
}