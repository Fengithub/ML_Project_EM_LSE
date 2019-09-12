library(reshape)
library(survival)
library(muhaz)
source('Seojin/myestimates.R')
source("Seojin/myerror_general.R")

my.linearerror.sub <- function(i, data, cvidx, modeltype, type, myerror = 1, mymax = Inf) {
  TRall <- data[c(cvidx$idx[cvidx$group != i]),]
  TR <- TRall[which(TRall$outcome == 1),]
  TE <- data[c(cvidx$idx[cvidx$group == i]),]
  if (type == "gaussian") {
    train.fit.all <- lm(log(time) ~ . , data = TR[,3:ncol(data)]) 
  } else if (type == "gamma_inverse") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "inverse"))
  } else if (type == "gamma_log") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "log"))
  } else if (type == "gamma_identity") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "identity"))
  } else {
    stop("Put your model type : gaussian, gamma_inverse, gamma_log, or gamma_identity?")
  }
  
  if (modeltype == "NoVar") {
    train.fit <- NULL
  } else if (modeltype == "AIC") {
    train.fit <- step(train.fit.all)
  } else if (modeltype == "BIC") {
    train.fit <- step(train.fit.all, k = nrow(TR))
  } else if (modeltype == "Full") {
    train.fit <- train.fit.all
  } else {
    stop("Put your model type : AIC or BIC or Full or NoVar?")
  }
  
if (modeltype != "NoVar") {
  if (type == "gaussian") {
    tnew <- exp(predict(train.fit, newdata = TE))
  } else if (type == "gamma_inverse") {
    tnew <- 1/(predict(train.fit, newdata = TE))
  } else if (type == "gamma_log") {
    tnew <- exp(predict(train.fit, newdata = TE))
  } else if (type == "gamma_identity") {
    tnew <- predict(train.fit, newdata = TE)
  } else {
    stop("Put your model type : gaussian, gamma_inverse, gamma_log, or gamma_identity?")
  }
  tnew <- sapply(tnew, min, mymax)
} else if (modeltype == "NoVar") {
  tnew <- min(sum(TE$time)/sum(TE$outcome), mymax)
}

  if (sum(TE$outcome)==0) {
  	err.nonrecur = sum(sapply((TE$time-tnew) * (1-TE$outcome), max, 0))/sum(1-TE$outcome) # nonrecur
  	err.recur = NA
  	err.all = err.nonrecur
  } else if (sum(1-TE$outcome)==0) {
  	err.nonrecur = NA
	err.recur.temp = (tnew-TE$time) * TE$outcome
	err.recur.temp[which(err.recur.temp >= -1 & err.recur.temp <= 0)] = 0
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

my.linearerror.sub.final <- function(i, data, cvidx, modeltype, type, myerror = 1, mymax = Inf) {
  TRall <- data[c(cvidx$idx[cvidx$group != i]),]
  TR <- TRall[which(TRall$outcome == 1),]
  TE <- data[c(cvidx$idx[cvidx$group == i]),]
  if (type == "gaussian") {
    train.fit.all <- lm(log(time) ~ . , data = TR[,3:ncol(data)]) 
  } else if (type == "gamma_inverse") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "inverse"))
  } else if (type == "gamma_log") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "log"))
  } else if (type == "gamma_identity") {
    train.fit.all <- glm(time ~ ., data = TR[,3:ncol(data)], family = Gamma(link = "identity"))
  } else {
    stop("Put your model type : gaussian, gamma_inverse, gamma_log, or gamma_identity?")
  }
  
  if (modeltype == "NoVar") {
    train.fit <- NULL
  } else if (modeltype == "AIC") {
    train.fit <- step(train.fit.all)
  } else if (modeltype == "BIC") {
    train.fit <- step(train.fit.all, k = nrow(TR))
  } else if (modeltype == "Full") {
    train.fit <- train.fit.all
  } else {
    stop("Put your model type : AIC or BIC or Full or NoVar?")
  }
  
if (modeltype != "NoVar") {
  if (type == "gaussian") {
    tnew <- exp(predict(train.fit, newdata = TE))
  } else if (type == "gamma_inverse") {
    tnew <- 1/(predict(train.fit, newdata = TE))
  } else if (type == "gamma_log") {
    tnew <- exp(predict(train.fit, newdata = TE))
  } else if (type == "gamma_identity") {
    tnew <- predict(train.fit, newdata = TE)
  } else {
    stop("Put your model type : gaussian, gamma_inverse, gamma_log, or gamma_identity?")
  }
  tnew <- sapply(tnew, min, mymax)
} else if (modeltype == "NoVar") {
  tnew <- min(sum(TE$time)/sum(TE$outcome), mymax)
}

  err.res <- myerrorfinal(predicted = tnew, outcome = TE$outcome, time = TE$time)
  return(err.res)
  #res <- data.frame(predicted = tnew, outcome = TE$outcome, time = TE$time, group = i, ID = TE$ID)
  #return(res)
}

my.linearerror <- function(data, cvidx, modeltype, type, mymax = Inf) {
  res <- t(sapply(1:10, my.linearerror.sub, data = datp, cvidx = myidx, modeltype = modeltype, type = type, mymax = mymax))
  colnames(res) <- c("All", "NonRecur", "Recur")
  return(res)
  }

my.linearpredicted <- function(data, cvidx, modeltype, type) {
  res <- lapply(1:10, my.linearerror.sub, data = datp, cvidx = myidx, modeltype = modeltype, type = type, myerror = 0)
  myres <- merge_recurse(res)
  return(myres)
  }
  
my.linearerror.final <- function(data, cvidx, modeltype, type, mymax = Inf) {
  res <- t(sapply(1:10, my.linearerror.sub.final, data = datp, cvidx = myidx, modeltype = modeltype, type = type, mymax = mymax))
  colnames(res) <- c("All", "NonRecur", "Recur")
  return(res)
  }


