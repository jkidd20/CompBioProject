### John Work File
source('~/RProjects/CompBioProject/Function File.R')
library(curatedBladderData)

data(GSE1827_eset)
a1 = GSE1827_eset

#geneCV(expData, label, pCut = .000001, nFold = 5, fitMethod = "rf", 
#           impMethod = "knnImpute", randomFlag = FALSE, nSelect = NA)

wholeTime = proc.time()
set.seed(2017)
# p value cut offs
tTime = Sys.time()
ss.p1.rf = geneCV(a1, pData(a1)$summarystage)
Sys.time() - tTime
# tTime = Sys.time()
# ss.p1.ada = geneCV(a1, pData(a1)$summarystage, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.p1.glmb = geneCV(a1, pData(a1)$summarystage, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.p1.adab = geneCV(a1, pData(a1)$summarystage, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.p2.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001)
Sys.time() - tTime
#tTime = Sys.time()
#ss.p2.ada = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001, fitMethod = "ada")
#Sys.time() - tTime
tTime = Sys.time()
ss.p2.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.p2.adab = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.p3.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001)
Sys.time() - tTime
#tTime = Sys.time()
#ss.p3.ada = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001, fitMethod = "ada")
#Sys.time() - tTime
tTime = Sys.time()
ss.p3.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.p3.adab = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.p4.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.001)
Sys.time() - tTime
#tTime = Sys.time()
#ss.p4.ada = geneCV(a1, pData(a1)$summarystage, pCut = 0.001, fitMethod = "ada")
#Sys.time() - tTime
tTime = Sys.time()
ss.p4.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.001, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.p4.adab = geneCV(a1, pData(a1)$summarystage, pCut = 0.001, fitMethod = "adaboost")
# Sys.time() - tTime

### Most Significant
tTime = Sys.time() 
ss.ns1.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 50)
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns1.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 50, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.ns1.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 50, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns1.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 50, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.ns2.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 100)
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns2.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 100, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.ns2.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 100, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns2.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 100, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.ns3.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 200)
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns3.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 200, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.ns3.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 200, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.ns3.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 200, fitMethod = "adaboost")
# Sys.time() - tTime

# Randomly Selected
tTime = Sys.time()
ss.rs1.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 50, randomFlag = TRUE)
Sys.time() - tTime
# tTime = Sys.time()
# ss.rs1.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 50, 
#                     randomFlag = TRUE, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.rs1.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 50,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.rs1.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 50,
#                      randomFlag = TRUE, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.rs2.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 100)
Sys.time() - tTime
#tTime = Sys.time()
#ss.rs2.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 100,
#                    randomFlag = TRUE, fitMethod = "ada")
#Sys.time() - tTime
tTime = Sys.time()
ss.rs2.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 100,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime
# tTime = Sys.time()
# ss.rs2.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 100,
#                      randomFlag = TRUE, fitMethod = "adaboost")
# Sys.time() - tTime

tTime = Sys.time()
ss.rs3.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 200, randomFlag = TRUE)
Sys.time() - tTime
# tTime = Sys.time()
# ss.rs3.ada = geneCV(a1, pData(a1)$summarystage, nSelect = 200,
#                     randomFlag = TRUE, fitMethod = "ada")
# Sys.time() - tTime
tTime = Sys.time()
ss.rs3.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 200,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime
#tTime = Sys.time()
#ss.rs3.adab = geneCV(a1, pData(a1)$summarystage, nSelect = 200,
#                     randomFlag = TRUE, fitMethod = "adaboost")
#Sys.time() - tTime

fullTime = wholeTime = proc.time()

#### heatmaps
