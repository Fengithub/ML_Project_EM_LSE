library(survival)
library(survAUC)
my.lp.auc.sub <- function(i, data, cvidx, modeltype, mytype, type) {
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
  lp <- predict(train.fit)
  lpnew <- predict(train.fit, newdata = TE)
  Surv.rsp <- Surv(TR$time, TR$outcome)
  Surv.rsp.new <- Surv(TE$time, TE$outcome)
  times <- seq(min(data$time), max(data$time)-1, length = 20)
  AUC_sh <- AUC.sh(Surv.rsp, Surv.rsp.new, lp, lpnew, times) 
  AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
  AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
  auclp.res <- c(AUC_sh$iauc, AUC_hc$iauc, AUC_Uno$iauc)
  return(auclp.res)
}

my.lp.auc <- function(data, cvidx, modeltype, type) {
  res <- t(sapply(1:10, my.lp.auc.sub, data = datp, cvidx = myidx, modeltype = modeltype, type = type))
  colnames(res) <- c("sh", "hc", "Uno")
  return(res)
  }
