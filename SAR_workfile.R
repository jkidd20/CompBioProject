## Sarah's Exploratory Work File

## Load appropriate libraries and dataset
library(BiocInstaller)
library(curatedBladderData)
data(GSE1827_eset)
a1 = GSE1827_eset
set.seed(2017)

##############################################################################
########### Changing the training set size ###################################
#############           Also adding smaller p-value cutoff           #########

## Finding differentially expressed genes on training sample of subset patients
ts20 = diffExpTable(20, a1)
ts30 = diffExpTable(30, a1)
ts40 = diffExpTable(40, a1)

## Prediction data sets for training size of 20
toPred = a1[, ts20$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage

## Gene sets to use for training size of 20
bigSet.g = rownames(ts20$geneTable)[ts20$geneTable$adj.P.Val < 0.001 & !is.na(ts20$geneTable$adj.P.Val)]
smallSet.g = rownames(ts20$geneTable)[ts20$geneTable$adj.P.Val < 0.0001 & !is.na(ts20$geneTable$adj.P.Val)]
smallerSet.g = rownames(ts20$geneTable)[ts20$geneTable$adj.P.Val < 0.00001 & !is.na(ts20$geneTable$adj.P.Val)]

## Impute missing values for training size of 20
imputedV = preProcess(toPred.expr.t, method = "knnImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

## Prediction testing for training size of 20
bigError.ts20.knn = cvError(bigSet.g, toPred.stage, toPred.expr.imp)
smallError.ts20.knn = cvError(smallSet.g, toPred.stage, toPred.expr.imp)
smallerError.ts20.knn = cvError(smallerSet.g, toPred.stage, toPred.expr.imp)

###############################################################################

## Prediction data sets for training size of 40
toPred = a1[, ts40$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage

## Gene sets to use for training size of 40
bigSet.g = rownames(ts40$geneTable)[ts40$geneTable$adj.P.Val < 0.001 & !is.na(ts40$geneTable$adj.P.Val)]
smallSet.g = rownames(ts40$geneTable)[ts40$geneTable$adj.P.Val < 0.0001 & !is.na(ts40$geneTable$adj.P.Val)]
smallerSet.g = rownames(ts40$geneTable)[ts40$geneTable$adj.P.Val < 0.00001 & !is.na(ts40$geneTable$adj.P.Val)]

## Impute missing values for training size of 40
imputedV = preProcess(toPred.expr.t, method = "knnImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

## Prediction testing for training size of 40
bigError.ts40.knn = cvError(bigSet.g, toPred.stage, toPred.expr.imp)
smallError.ts40.knn = cvError(smallSet.g, toPred.stage, toPred.expr.imp)
smallerError.ts40.knn = cvError(smallerSet.g, toPred.stage, toPred.expr.imp)

###############################################################################
############# Changing the imputation method; train size = 30 #################
#############           Also adding smaller p-value cutoff           ##########

## Prediction data sets for training size of 30
toPred = a1[, ts30$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage

## Gene sets to use for training size of 30
bigSet.g = rownames(ts30$geneTable)[ts30$geneTable$adj.P.Val < 0.001 & !is.na(ts30$geneTable$adj.P.Val)]
smallSet.g = rownames(ts30$geneTable)[ts30$geneTable$adj.P.Val < 0.0001 & !is.na(ts30$geneTable$adj.P.Val)]
smallerSet.g = rownames(ts30$geneTable)[ts30$geneTable$adj.P.Val < 0.00001 & !is.na(ts30$geneTable$adj.P.Val)]

## Impute missing values for training size of 30 - with MEDIAN method
imputedV = preProcess(toPred.expr.t, method = "medianImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

## Prediction testing for training size of 30
bigError.ts30.med = cvError(bigSet.g, toPred.stage, toPred.expr.imp)
smallError.ts30.med = cvError(smallSet.g, toPred.stage, toPred.expr.imp)
smallerError.ts30.med = cvError(smallerSet.g, toPred.stage, toPred.expr.imp)

###############################################################################
############### Using KNN imputation method; train size = 30 ##################
#######   Random choice of 50 differentially expressed genes       ############

## Prediction data sets for training size of 30
toPred = a1[, ts30$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage

## Randomly choose 50 genes for training size of 30
notNASet = ts30$geneTable[!is.na(ts30$geneTable$adj.P.Val), ]
rand50Set.g = rownames(notNASet[sample(nrow(notNASet), size=50), ])

## Impute missing values for training size of 30 - with MEDIAN method
imputedV = preProcess(toPred.expr.t, method = "medianImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

## Prediction testing for training size of 30
rand50Error.ts30.med = cvError(rand50Set.g, toPred.stage, toPred.expr.imp)


