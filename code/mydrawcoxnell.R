my.draw.coxsnell <- function(model, data) {
  source("ggsurv.R")
  library(survival)
  library(ggplot2)
  dfree = data$outcome
  mres = resid(model, type="martingale")
  csres = data$outcome-mres
  r.surv = survfit(Surv(csres,dfree)~1, type="fleming-harrington")
  ggplot(data = data.frame(y = r.surv$time, x = -log(r.surv$surv)), 
         mapping = aes(x = y, y = x)) +
    geom_step(direction = "hv") + 
    ylab("Estimated Cum Hazards") + xlab("Cox-Snell residual") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2)
}
