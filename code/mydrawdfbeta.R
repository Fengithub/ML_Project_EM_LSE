my.draw.dfbeta <- function(model) {
  library(survival)
  library(ggplot2)
  library(reshape)
  my.dfbeta <- residuals(model, type="dfbetas")
  colnames(my.dfbeta) = names(model$assign)
  my.dfbeta <- as.data.frame(my.dfbeta)
  my.dfbeta$idx <- 1:nrow(my.dfbeta)
  my.resid.dfbeta.melt <- melt(my.dfbeta, id.vars = c("idx"))
  ggplot(data = my.resid.dfbeta.melt, mapping = aes(x = idx, y = value)) +
    layer(stat = "identity", geom = "point", position = "identity") +
    facet_wrap( ~ variable, scales = "free") +
    scale_y_continuous(name = "DFBETA")
}
