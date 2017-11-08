# Playing Around
# To install and load data - Curated Bladder
library(curatedBladderData)

## Note on ages in GSE1827 - there is a patient that is supposedly 113 years old... Possibly remove?
data(GSE1827_eset)
a1 = GSE1827_eset

#### Finding differentially expressed genes on training sample of subset patients
set.seed(2017)
library(limma)

diffExpTable = function(trainSize = 30, expData){

  eSet = exprs(expData)
  
  invNum = which(pData(expData)$summarystage == "invasive")
  supNum = which(pData(expData)$summarystage == "superficial")
  
  trainSamp = c(sample(invNum, trainSize / 2), sample(supNum, trainSize / 2))
  predSamp = seq_len(ncol(expData))[!(seq_len(ncol(expData)) %in% trainSamp)]
  
  toTrain = expData[, trainSamp]
  
  tumType = pData(toTrain)$summarystage
  design <- model.matrix(~ tumType)
  train.eSet = exprs(toTrain)
  
  fit <- lmFit(train.eSet, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, n=nrow(eSet), sort.by="none")
  
  return(list(predSamp = predSamp, geneTable = tab))
}

try1 = diffExpTable(30, a1)

### Prediction data sets
toPred = a1[, try1$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage


### Gene sets to use
bigSet.g = rownames(try1$geneTable)[try1$geneTable$adj.P.Val < 0.001 & !is.na(try1$geneTable$adj.P.Val)]
smallSet.g = rownames(try1$geneTable)[try1$geneTable$adj.P.Val < 0.0001 & !is.na(try1$geneTable$adj.P.Val)]


# impute missing values
library(caret)
library(e1071)
imputedV = preProcess(toPred.expr.t, method = "knnImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

##### Prediction testing

cvError = function(geneSet, label, expSet, nFold = 10, fitMethod = "rf"){
  toSample = matrix(1:nFold, 1, length(label))
  cvNum = sample(toSample)
  numWrong = rep(0, nFold)
  for(i in 1:nFold){
      # set up Train set
      trainDF = data.frame(label[cvNum != i], t(expSet)[cvNum != i, rownames(expSet) %in% geneSet] )
      colnames(trainDF)[1] = "label"
       # train data
      fit = train(label ~ ., data = trainDF, method = fitMethod)
      
      # set up Test Set
      testDF = data.frame(label[cvNum == i], t(expSet)[cvNum == i, rownames(expSet) %in% geneSet] )
      colnames(testDF)[1] = "label"
      
      predValues = predict(fit, testDF)
      numWrong[i] = sum(as.character(predValues) != as.character(testDF$label))
  }
  mean(numWrong / table(cvNum))
}

bigError = cvError(bigSet.g, toPred.stage, toPred.expr.imp)
smallError = cvError(smallSet.g, toPred.stage, toPred.expr.imp)
