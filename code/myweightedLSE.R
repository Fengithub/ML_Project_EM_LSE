#source("myweightedLSE.R")

myweightedLSE <- function(time, outcome, train, hazard, survival, y = 1:130, type = "gaussian", modeltype = "Full") {
	idx = which(outcome == 1)	# index for recur data
	ttr = time[idx]
	dfs = time[-idx]
	train.ttr = train[idx,]
	train.dfs = train[-idx,]
	my.survival = survival[match(dfs, survival[,1]), 2]
	my.weight.ttr = rep(1, length(idx))
	my.weight.dfs0 = lapply(dfs, f <- function(x, y) return(x*y), y)
	my.weight.dfs0 = unlist(my.weight.dfs0)
	my.weight.dfs = hazard[match(my.weight.dfs0, hazard[,1]),2]
	my.weight.dfs[is.na(my.weight.dfs)] = 0
	newweight = c(my.weight.ttr, my.weight.dfs)
	TR = rbind(train.ttr, apply(train.dfs, 2, rep, each = length(y)))
	if (type == "gaussian") {
		newtime.dfs = rep(log(dfs), each = length(y)) + rep(log(y), length(dfs))
		newtime = c(log(ttr), newtime.dfs)
		my.fit = lm(newtime ~ as.matrix(TR), weights = newweight)
	} else if (type == "gamma_log") {
		newtime.dfs = rep(dfs, each = length(y)) + rep(y, length(dfs))
		newtime = c(ttr, newtime.dfs)
		my.fit = glm(newtime ~ TR, weights = newweight, family = Gamma(link = "log"))
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
  	train.fit <- step(my.fit, k = nrow(train))
		coefnew = names(train.fit$coefficients)
 	} else if (modeltype == "Full") {
  		train.fit <- my.fit
		coefnew = names(train.fit$coefficients)
 	} else {
		stop("Put your model type : AIC or BIC or Full or NoVar?")
 	}

	res = list(beta = train.fit$coefficients, coefficients = coefnew)
	return(res)
}