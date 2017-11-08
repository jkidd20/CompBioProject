######## File holding Functions
library(limma)
library(caret)
library(e1071)

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