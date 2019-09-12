my.draw.deviance <- function(model, data) {
  library(survival)
  library(reshape)
  library(ggplot2)
  my.res.dev = residuals(model, type = "deviance") # weird
  my.resid.dev <- data.frame(my.res.dev, row.names = NULL)
  my.resid.dev$LinearPredictor <- predict(model)
  my.resid.dev <- cbind(my.resid.dev, data[,names(model$assign)])
  my.resid.dev.melt <- melt(my.resid.dev, id.vars = c("my.res.dev"))
  ggplot(data = my.resid.dev.melt, mapping = aes(x = value, y = my.res.dev)) +
    layer(stat = "identity", geom = "point", position = "identity") +
    facet_wrap( ~ variable, scales = "free") +
    stat_smooth() + 
    scale_y_continuous(name = "Deviance residuals")
}
