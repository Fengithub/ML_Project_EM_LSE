#model=survfit(coxph.fit.subAIC,type='aalen')
library(ggplot2)
source("codes/multiplot.R")
mymodelchecking <- function(model) {
  my.para.asse <- data.frame(time = model$time, Exp =-log(model$surv),
                             Wei = log(-log(model$surv)),
                             LogLogistic = -log(model$surv/(1-model$surv)),
                             LogNormal = qnorm(1-model$surv))
  fit.exp <- lm(Exp ~ time, data = my.para.asse)
  p1 <- ggplot(data = my.para.asse, 
    mapping = aes(x = time, y = Exp)) +
    geom_step(direction = "hv") + 
    ylab("Estimated Cumulative Hazard Rates") + xlab("Time") +
    ggtitle("Exponential Model Check") + 
    geom_abline(slope = fit.exp$coefficients[2], intercept = fit.exp$coefficients[1], color = "red", linetype = 2)
  fit.wei <- lm(Wei ~ log(time), data = my.para.asse)
  p2 <- ggplot(data = my.para.asse, 
               mapping = aes(x = log(time), y = Wei)) +
    geom_step(direction = "hv") + 
    ylab("Log Estimated Cumulative Hazard Rates") + xlab("Log of Time") +
    ggtitle("Weibull Model Check") + 
    geom_abline(slope = fit.wei$coefficients[2], 
                intercept = fit.wei$coefficients[1], 
                color = "red", linetype = 2)
  fit.ll <- lm(LogLogistic ~ log(time), data = my.para.asse)
  p3 <- ggplot(data = my.para.asse, 
               mapping = aes(x = log(time), y = LogLogistic)) +
    geom_step(direction = "hv") + 
    ylab("Log of Odds") + xlab("Log of Time") +
    ggtitle("Log-Logistic Model Check") + 
    geom_abline(slope = fit.ll$coefficients[2], 
                intercept = fit.ll$coefficients[1], 
                color = "red", linetype = 2)
  fit.ln <- lm(LogNormal ~ log(time), data = my.para.asse)
  p4 <- ggplot(data = my.para.asse, 
               mapping = aes(x = log(time), y = LogNormal)) +
    geom_step(direction = "hv") + 
    ylab("Normal Quantile") + xlab("Log of Time") +
    ggtitle("Log-Normal Model Check") + 
    geom_abline(slope = fit.ln$coefficients[2], 
                intercept = fit.ln$coefficients[1], 
                color = "red", linetype = 2)
  r2 <- c(summary(fit.exp)$r.squared,
          summary(fit.wei)$r.squared,
          summary(fit.ll)$r.squared,
          summary(fit.ln)$r.squared)
  r2 <- as.data.frame(t(r2))
  colnames(r2) <- c("Exponential", "Weibull", "Log-Logistic", "Log-Normal")
  multiplot(p1, p2, p3, p4, cols = 2)
  return(r2)
}