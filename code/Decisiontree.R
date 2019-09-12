setwd('/Users/fenpei/Box Sync/Project/Data')
alldata<-read.table('wpbc.newdata.txt', header = T)
recurdata <- alldata[which(alldata$outcome == 1),]
nonrecurdata <- alldata[which(alldata$outcome == 0),]


set.seed(1234)
# sample train data set and test data set
ind<-sample(2,nrow(alldata),replace=TRUE,prob=c(0.7,0.3))
trainData<-alldata[ind==1,]
testData<-alldata[ind==2,]
# construct a tree
library(party)
library(survival)
my.surv <- Surv(alldata$time, alldata$outcome)
my.ctrl <-  ctree_control(teststat = c("quad", "max"),
  testtype = "MonteCarlo",
  mincriterion = 0.1, minsplit = 20, minbucket = 7,
  stump = FALSE, nresample = 9999, maxsurrogate = 0,
  mtry = 0, savesplitstats = TRUE, maxdepth = 0)

Recur_ctree<-ctree(my.surv ~ ., data=alldata[,-c(1:3)])
Recur_ctree<-ctree(my.surv ~ ., data=alldata[,-c(1:3)], controls = my.ctrl)
plot(Recur_ctree)



#make prediction
table(predict(Recur_ctree),trainData$time)
#print results
print(Recur_ctree)
plot(Recur_ctree)
plot(Recur_ctree,type='simple')
#make prediction on test data
testPred<-predict(Recur_ctree,newdata=testData)
table(testPred,testData$time)



