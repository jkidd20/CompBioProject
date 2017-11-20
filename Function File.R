######## File holding Functions
library(limma)
library(caret)
library(e1071)
library(pheatmap)
library(RColorBrewer)

geneCV = function(expData, label, pCut = .000001, nFold = 5, fitMethod = "rf", impMethod = "medianImpute", 
                  randomFlag = FALSE, nSelect = NA){
    eSet.t = exprs(expData)
    label = as.factor(label)
    labLev = levels(label)
    g1Num = which(label == labLev[1])
    g2Num = which(label == labLev[2])
    
    ### Imputation
    imputedV = preProcess(t(eSet.t), method = impMethod)
    eSet = t(predict(imputedV, t(eSet.t)))
    
    ### Cross Validation
    toSample1 = matrix(1:nFold, 1, length(g1Num))
    cvNum1 = sample(toSample1)
    toSample2 = matrix(1:nFold, 1, length(g2Num))
    cvNum2 = sample(toSample2)
    
    numWrong = rep(0, nFold)
    genes = list()
    
    for(i in 1:nFold){
      # set up Train and Test set
      trainLabels = c(rep(labLev[1], sum(cvNum1 != i)), rep(labLev[2], sum(cvNum2 != i)))
      testLabels = c(rep(labLev[1], sum(cvNum1 == i)), rep(labLev[2], sum(cvNum2 == i)))
      
      trainExprs = rbind(t(eSet)[g1Num[cvNum1 != i], ], t(eSet)[g2Num[cvNum2 != i], ])
      testExprs = rbind(t(eSet)[g1Num[cvNum1 == i], ], t(eSet)[g2Num[cvNum2 == i], ])
      
      
      trainDF = data.frame(trainLabels, trainExprs)
      testDF = data.frame(testLabels, testExprs)
      
      colnames(trainDF)[1] = colnames(testDF)[1] = "label"
      
      # train data
      if(!randomFlag){
        design <- model.matrix(~ trainDF[, 1])
        fit <- lmFit(t(trainExprs), design)
        fit <- eBayes(fit)
        tab <- topTable(fit, n=nrow(eSet), sort.by="none")
        if(is.na(nSelect)){
          genes[[i]] = gsub('-', '.', rownames(tab)[tab$adj.P.Val < pCut])
        }
        else{
          genes[[i]] = gsub('-', '.', rownames(tab[order(tab$adj.P.Val),])[1:nSelect])
        }
      }
      else{
        genes[[i]] = gsub('-', '.', sample(rownames(expData), nSelect))
      }
      
      smallTrain = cbind(trainDF[, 1], trainDF[, colnames(trainDF) %in% genes[[i]]])
      colnames(smallTrain)[1] = "label"
      
      mFit = train(label ~ ., data = smallTrain, method = fitMethod)
      
      # Predict on Test
      
      predValues = predict(mFit, testDF)
      numWrong[i] = sum(as.character(predValues) != as.character(testDF$label))
    }
    return(list(genes = genes, cvError = mean(numWrong / table(c(cvNum1, cvNum2))), 
                impMethod = impMethod, label = label))
}


#### Function to create heatmaps from first function
myMap = function(cvRes, expData, titleIn = "", occur = 1){
  plotGenes = names(table(unlist(cvRes$genes)))[table(unlist(cvRes$genes)) >= occur]

  eSet.t = exprs(expData)
  label = cvRes$label
  labLev = levels(label)
  
  ### Imputation
  imputedV = preProcess(t(eSet.t), method = cvRes$impMethod)
  eSet = t(predict(imputedV, t(eSet.t)))
  
  dists = dist(t(eSet))
  cols <- colorRampPalette(rev(brewer.pal(9,"RdBu")))(255)
  distmat <- as.matrix(dists)
  df <- data.frame(condition=label,
                   row.names=colnames(distmat))
  
  eSet.plot = eSet[gsub('-', '.', rownames(expData)) %in% plotGenes, ]

  pheatmap(eSet.plot, color=cols,
           annotation_col=df,
           show_rownames=FALSE, show_colnames=FALSE)
  }

# ##### Old Functions
# diffExpTable = function(trainSize = 30, expData){
#   
#   eSet = exprs(expData)
#   
#   invNum = which(pData(expData)$summarystage == "invasive")
#   supNum = which(pData(expData)$summarystage == "superficial")
#   
#   trainSamp = c(sample(invNum, trainSize / 2), sample(supNum, trainSize / 2))
#   predSamp = seq_len(ncol(expData))[!(seq_len(ncol(expData)) %in% trainSamp)]
#   
#   toTrain = expData[, trainSamp]
#   
#   tumType = pData(toTrain)$summarystage
#   design <- model.matrix(~ tumType)
#   train.eSet = exprs(toTrain)
#   
#   fit <- lmFit(train.eSet, design)
#   fit <- eBayes(fit)
#   tab <- topTable(fit, n=nrow(eSet), sort.by="none")
#   
#   return(list(predSamp = predSamp, geneTable = tab))
# }
# 
# ##### Prediction testing
# 
# cvError = function(geneSet, label, expSet, nFold = 10, fitMethod = "rf"){
#   toSample = matrix(1:nFold, 1, length(label))
#   cvNum = sample(toSample)
#   numWrong = rep(0, nFold)
#   for(i in 1:nFold){
#     # set up Train set
#     trainDF = data.frame(label[cvNum != i], t(expSet)[cvNum != i, rownames(expSet) %in% geneSet] )
#     colnames(trainDF)[1] = "label"
#     # train data
#     fit = train(label ~ ., data = trainDF, method = fitMethod)
#     
#     # set up Test Set
#     testDF = data.frame(label[cvNum == i], t(expSet)[cvNum == i, rownames(expSet) %in% geneSet] )
#     colnames(testDF)[1] = "label"
#     
#     predValues = predict(fit, testDF)
#     numWrong[i] = sum(as.character(predValues) != as.character(testDF$label))
#   }
#   mean(numWrong / table(cvNum))
# }