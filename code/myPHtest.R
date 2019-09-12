my.PHtest <- function(model, data, threshold = 0.05, filename = 'full', plot = 'False') {
    library(survival)
    model.ptest <- cox.zph(model)		## proportionality test
  if (plot == 'True') {
    for (i in names(model$assign)) {
      myVarName = colnames(data)[i]
      #png(paste('SurvivalAnalysis_CoxPH_', filename, '_PropAssumpTest_', i, '.png', sep = ""))
      plot(model.ptest, var = i)
      #dev.off()
    }
  }	
  my.res <- as.data.frame(model.ptest$table)
  my.res$result <- "OK"
  my.res$result[model.ptest$table[,"p"] < threshold] <- "Not_OK"
  return(my.res)
}