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
tTime = Sys.time()
ss.p1.glmb = geneCV(a1, pData(a1)$summarystage, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.p2.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001)
Sys.time() - tTime
tTime = Sys.time()
ss.p2.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.00001, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.p3.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001)
Sys.time() - tTime
tTime = Sys.time()
ss.p3.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.0001, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.p4.rf = geneCV(a1, pData(a1)$summarystage, pCut = 0.001)
Sys.time() - tTime
tTime = Sys.time()
ss.p4.glmb = geneCV(a1, pData(a1)$summarystage, pCut = 0.001, fitMethod = "glmboost")
Sys.time() - tTime

### Most Significant
tTime = Sys.time() 
ss.ns1.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 50)
Sys.time() - tTime
tTime = Sys.time()
ss.ns1.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 50, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.ns2.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 100)
Sys.time() - tTime
tTime = Sys.time()
ss.ns2.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 100, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.ns3.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 200)
Sys.time() - tTime
tTime = Sys.time()
ss.ns3.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 200, fitMethod = "glmboost")
Sys.time() - tTime

# Randomly Selected
tTime = Sys.time()
ss.rs1.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 50, randomFlag = TRUE)
Sys.time() - tTime
tTime = Sys.time()
ss.rs1.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 50,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.rs2.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 100)
Sys.time() - tTime
tTime = Sys.time()
ss.rs2.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 100,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime

tTime = Sys.time()
ss.rs3.rf = geneCV(a1, pData(a1)$summarystage, nSelect = 200, randomFlag = TRUE)
Sys.time() - tTime
tTime = Sys.time()
ss.rs3.glmb = geneCV(a1, pData(a1)$summarystage, nSelect = 200,
                     randomFlag = TRUE, fitMethod = "glmboost")
Sys.time() - tTime

fullTime = wholeTime = proc.time()

#### heatmaps
myMap(ss.p1.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.00001) - Occur Once", occur = 1)
myMap(ss.p1.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.00001) - Occurs Twice", occur = 2)
myMap(ss.p1.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.00001) - Occurs in All", occur = 5)

myMap(ss.p2.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.0001) - Occur Once", occur = 1)
myMap(ss.p2.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.0001) - Occurs Twice", occur = 2)
myMap(ss.p2.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.0001) - Occurs in All", occur = 5)

myMap(ss.p3.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.001) - Occur Once", occur = 1)
myMap(ss.p3.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.001) - Occurs Twice", occur = 2)
myMap(ss.p3.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.001) - Occurs in All", occur = 5)

myMap(ss.p4.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.01) - Occur Once", occur = 1)
myMap(ss.p4.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.01) - Occurs Twice", occur = 2)
myMap(ss.p4.rf, a1, titleIn = "Differentially Expressed Genes (p < 0.01) - Occurs in All", occur = 5)

myMap(ss.ns1.rf, a1, titleIn = "Most Significant Genes - Top 50 - Occur Once", occur = 1)
myMap(ss.ns1.rf, a1, titleIn = "Most Significant Genes - Top 50 - Occurs Twice", occur = 2)
myMap(ss.ns1.rf, a1, titleIn = "Most Significant Genes - Top 50 - Occurs in All", occur = 5)

myMap(ss.ns2.rf, a1, titleIn = "Most Significant Genes - Top 100 - Occur Once", occur = 1)
myMap(ss.ns2.rf, a1, titleIn = "Most Significant Genes - Top 100 - Occurs Twice", occur = 2)
myMap(ss.ns2.rf, a1, titleIn = "Most Significant Genes - Top 100 - Occurs in All", occur = 5)

myMap(ss.ns3.rf, a1, titleIn = "Most Significant Genes - Top 200 - Occur Once", occur = 1)
myMap(ss.ns3.rf, a1, titleIn = "Most Significant Genes - Top 200 - Occurs Twice", occur = 2)
myMap(ss.ns3.rf, a1, titleIn = "Most Significant Genes - Top 200 - Occurs in All", occur = 5)

myMap(ss.rs1.rf, a1, titleIn = "Randomly Selected Genes - 50 - Occur Once", occur = 1)
myMap(ss.rs1.rf, a1, titleIn = "Randomly Selected Genes - 50 - Occurs Twice", occur = 2)
myMap(ss.rs1.rf, a1, titleIn = "Randomly Selected Genes - 50 - Occurs in All", occur = 5)

myMap(ss.rs2.rf, a1, titleIn = "Randomly Selected Genes - 100 - Occur Once", occur = 1)
myMap(ss.rs2.rf, a1, titleIn = "Randomly Selected Genes - 100 - Occurs Twice", occur = 2)
myMap(ss.rs2.rf, a1, titleIn = "Randomly Selected Genes - 100 - Occurs in All", occur = 5)

myMap(ss.rs3.rf, a1, titleIn = "Randomly Selected Genes - 200 - Occur Once", occur = 1)
myMap(ss.rs3.rf, a1, titleIn = "Randomly Selected Genes - 200 - Occurs Twice", occur = 2)
myMap(ss.rs3.rf, a1, titleIn = "Randomly Selected Genes - 200 - Occurs in All", occur = 5)