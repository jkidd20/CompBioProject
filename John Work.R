### John Work File
source('~/RProjects/CompBioProject/Function File.R')
library(curatedBladderData)

data(GSE1827_eset)
a1 = GSE1827_eset

johnTry = diffExpTable(30, a1)

toPred = a1[, johnTry$predSamp]
toPred.expr.t = data.frame(exprs(toPred))
toPred.stage = pData(toPred)$summarystage


### Gene sets to use
gSet1 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 0.00001 & !is.na(johnTry$geneTable$adj.P.Val)]
gSet2 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 0.0001 & !is.na(johnTry$geneTable$adj.P.Val)]
gSet3 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 0.001 & !is.na(johnTry$geneTable$adj.P.Val)]
gSet4 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 0.01 & !is.na(johnTry$geneTable$adj.P.Val)]
gSet5 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 0.1 & !is.na(johnTry$geneTable$adj.P.Val)]
gSet6 = rownames(johnTry$geneTable)[johnTry$geneTable$adj.P.Val < 1.1 & !is.na(johnTry$geneTable$adj.P.Val)]


# impute missing values
library(caret)
library(e1071)
imputedV = preProcess(toPred.expr.t, method = "knnImpute")
toPred.expr.imp = predict(imputedV, toPred.expr.t)

#Random Forests
cvError(gSet1, toPred.stage, toPred.expr.imp)
cvError(gSet2, toPred.stage, toPred.expr.imp)
cvError(gSet3, toPred.stage, toPred.expr.imp)
cvError(gSet4, toPred.stage, toPred.expr.imp)
cvError(gSet5, toPred.stage, toPred.expr.imp)
cvError(gSet6, toPred.stage, toPred.expr.imp)

# logistic regression
cvError(gSet1, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
cvError(gSet2, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
cvError(gSet3, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
cvError(gSet4, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
cvError(gSet5, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
cvError(gSet6, toPred.stage, toPred.expr.imp, fitMethod = "glmboost")
