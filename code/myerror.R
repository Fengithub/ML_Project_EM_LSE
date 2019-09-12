library(survival)
library(survAUC)
my.error.sub0 <- function(i, data, cvidx, modeltype, mytype, type) {
  TR <- data[c(cvidx$idx[cvidx$group != i]),]
  TE <- data[c(cvidx$idx[cvidx$group == i]),]
  train.surv <- Surv(TR$time, TR$outcome)
  if (type == "cox") {
    train.fit.all <- coxph(train.surv ~ . , method="breslow", 
                           data = TR[,4:ncol(data)]) 
  } else if (type == "exponential") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "exponential")
  } else if (type == "gaussian") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "gaussian")
  } else if (type == "weibull") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "weibull")
  } else if (type == "lognormal") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "lognormal")
  } else if (type == "logistic") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "logistic")
  } else if (type == "loglogistic") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "loglogistic")
  } else {
    stop("Put your model type : cox, exponential, weibull, gaussian, lognormal, loglogstic, or loglogistic?")
  }
  
  if (modeltype == "AIC") {
    train.fit <- step(train.fit.all)
  } else if (modeltype == "BIC") {
    train.fit <- step(train.fit.all, k = nrow(TR))
  } else if (modeltype == "Full") {
    train.fit <- train.fit.all
  } else {
    stop("Put your model type : AIC or BIC or Full?")
  }
  tnew <- predict(train.fit, newdata = TE,  type="response")

  if (sum(TE$outcome)==0) {
  	err.nonrecur = sum(sapply((TE$time-tnew) * (1-TE$outcome), max, 0))/sum(1-TE$outcome) # nonrecur
  	err.recur = NA
  	err.all = err.nonrecur
  } else if (sum(1-TE$outcome)==0) {
  	err.nonrecur = NA
  	err.recur = sum(sapply((tnew-TE$time) * TE$outcome, max, 0))/sum(TE$outcome) # recur
  	err.all = err.nonrecur
  } else {
  	err.recur = sum(sapply((tnew-TE$time) * TE$outcome, max, 0))/sum(TE$outcome) # recur
	err.nonrecur = sum(sapply((TE$time-tnew) * (1-TE$outcome), max, 0))/sum(1-TE$outcome) # nonrecur
	err.all = (err.recur*sum(TE$outcome) + err.nonrecur*sum(1-TE$outcome))/nrow(TE)
  }
  err.res <- c(err.all, err.nonrecur, err.recur)
  return(err.res)
}

my.error.sub <- function(i, data, cvidx, modeltype, mytype, type, myerror = 1, mymax = Inf) {
  TR <- data[c(cvidx$idx[cvidx$group != i]),]
  TE <- data[c(cvidx$idx[cvidx$group == i]),]
  train.surv <- Surv(TR$time, TR$outcome)
  if (type == "cox") {
    train.fit.all <- coxph(train.surv ~ . , method="breslow", 
                           data = TR[,4:ncol(data)]) 
  } else if (type == "exponential") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "exponential")
  } else if (type == "gaussian") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "gaussian")
  } else if (type == "weibull") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "weibull")
  } else if (type == "lognormal") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "lognormal")
  } else if (type == "logistic") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "logistic")
  } else if (type == "loglogistic") {
    train.fit.all <- survreg(train.surv ~ ., data = TR[,4:ncol(data)], dist = "loglogistic")
  } else {
    stop("Put your model type : cox, exponential, weibull, gaussian, lognormal, loglogstic, or loglogistic?")
  }
  
  if (modeltype == "AIC") {
    train.fit <- step(train.fit.all)
  } else if (modeltype == "BIC") {
    train.fit <- step(train.fit.all, k = nrow(TR))
  } else if (modeltype == "Full") {
    train.fit <- train.fit.all
  } else {
    stop("Put your model type : AIC or BIC or Full?")
  }
  tnew <- predict(train.fit, newdata = TE,  type="response")
  tnew <- sapply(tnew, min, mymax)
  
  if (sum(TE$outcome)==0) {
  	err.nonrecur = sum(sapply((TE$time-tnew) * (1-TE$outcome), max, 0))/sum(1-TE$outcome) # nonrecur
  	err.recur = NA
  	err.all = err.nonrecur
  } else if (sum(1-TE$outcome)==0) {
  	err.nonrecur = NA
  	err.recur.temp = (tnew-TE$time) * TE$outcome
	err.recur.temp[which(err.recur.temp == -1)] = 0
  	err.recur = sum(abs(err.recur.temp))/sum(TE$outcome) # recur
  	err.all = err.nonrecur
  } else {
  	err.recur.temp = (tnew-TE$time) * TE$outcome
	err.recur.temp[which(err.recur.temp >= -1 & err.recur.temp <= 0)] = 0
  	err.recur = sum(abs(err.recur.temp))/sum(TE$outcome) # recur
  	err.nonrecur = sum(sapply((TE$time-tnew) * (1-TE$outcome), max, 0))/sum(1-TE$outcome) # nonrecur
	err.all = (err.recur*sum(TE$outcome) + err.nonrecur*sum(1-TE$outcome))/nrow(TE)
  }
  if (myerror == 1) {
  	err.res <- c(err.all, err.nonrecur, err.recur)
  	} else {
  	err.res <- data.frame(predicted = tnew, time = TE$time, outcome = TE$outcome)
  	}
  return(err.res)
}

my.error0 <- function(data, cvidx, modeltype, type) {
  res <- t(sapply(1:10, my.error.sub0, data = datp, cvidx = myidx, modeltype = modeltype, type = type))
  colnames(res) <- c("All", "NonRecur", "Recur")
  return(res)
  }

my.error <- function(data, cvidx, modeltype, type, mymax = Inf) {
  res <- t(sapply(1:10, my.error.sub, data = datp, cvidx = myidx, modeltype = modeltype, type = type, mymax = mymax))
  colnames(res) <- c("All", "NonRecur", "Recur")
  return(res)
  }
  
my.linearpredicted <- function(data, cvidx, modeltype, type) {
  res <- lapply(1:10, my.error.sub, data = datp, cvidx = myidx, modeltype = modeltype, type = type, myerror = 0)
  myres <- merge_recurse(res)
  return(myres)
  }

