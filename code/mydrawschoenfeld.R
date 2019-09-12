my.draw.schoenfeld <- function(model) {
  library(survival)
  library(ggplot2)
  library(reshape)
  my.res.sch = residuals(model, type = "schoenfeld")
  my.survtime <- as.numeric(rownames(my.res.sch))
  my.resid.schoenfeld <- data.frame(my.res.sch, row.names = NULL)
  my.resid.schoenfeld$survtime <- my.survtime
  my.resid.schoenfeld.melt <- melt(my.resid.schoenfeld, id.vars = c("survtime"))
  ggplot(data = my.resid.schoenfeld.melt, mapping = aes(x = survtime, y = value)) +
    layer(stat = "identity", geom = "point", position = "identity") +
    facet_wrap( ~ variable, scales = "free") +
    stat_smooth() + 
    scale_y_continuous(name = "Schoenfeld residuals")
}
