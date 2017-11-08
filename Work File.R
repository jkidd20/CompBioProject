# Playing Around
# To install and load data - Curated Bladder
library(BiocInstaller)
biocLite("curatedBladderData")
library(curatedBladderData)

browseVignettes("curatedBladderData")
data(package="curatedBladderData")

### Which dataset to use? 

##########################################
### GSE31189 has the most phenotypic data

# Link to the reference paper
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3537330/pdf/nihms415409.pdf

###################
# Paper doesn't have much that is terribly interesting,
# but data was used for biomarkers, could be really good 
# for clustering, meat maps, PCA?
# Will definitely need to do some adjustments for multiple testing

#### Possible approach
# Determine which genes are differentially expressed
# Use PCA and clustering in EDA for this approach
# Run tests, correction
# plot heat map of DE genes

############
# Paper focuses on cancer and healthy subjects, but we could have some fun and look
# at some differences by the other variables, such as gender and what not as a step
# in determining which are significant between cancers.

# Load data
# Cancer/No Cancer - GSE31189

## Note on ages in GSE1827 - there is a patient that is supposedly 113 years old... Possibly remove?
data(GSE1827_eset)
a1 = GSE1827_eset

####### Suggested to use limma (lmFit, etc) to fit model, find prediction
### Mike will give us code
set.seed(2017)
library(limma)

eSet = exprs(a1)

invNum = which(pData(a1)$summarystage == "invasive")
supNum = which(pData(a1)$summarystage == "superficial")

trainSamp = c(sample(invNum, 15), sample(supNum, 15))
predSamp = seq_len(ncol(a1))[!(seq_len(ncol(a1)) %in% trainSamp)]

toTrain = a1[, trainSamp]

tumType = pData(toTrain)$summarystage
design <- model.matrix(~ tumType)
train.eSet = exprs(toTrain)

fit <- lmFit(train.eSet, design)
fit <- eBayes(fit)
tab <- topTable(fit, n=nrow(eSet), sort.by="none")


bigSet.g = rownames(tab)[tab$adj.P.Val < 0.001 & !is.na(tab$adj.P.Val)]
smallSet.g = rownames(tab)[tab$adj.P.Val < 0.0001 & !is.na(tab$adj.P.Val)]


toPred = a1[, predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage
# impute missing values
library(caret)
library(e1071)
imputedV = preProcess(toPred.expr.t, method = "knnImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

toPred.bigexp = exprs(toPred)[rownames(toPred) %in% bigSet.g, ] 
toPred.smallexp = exprs(toPred)[rownames(toPred) %in% smallSet.g, ] 

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
