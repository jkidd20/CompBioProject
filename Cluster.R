# Recreate Clustering Diagram
library(curatedBladderData)
library(dendextend)
library(RColorBrewer)
library(matrixStats)
data("GSE1827_eset")
a1 = GSE1827_eset

eSet = exprs(a1)
numMis = rowSums(is.na(eSet))
eSet2 = eSet[numMis <= 2, ]
dists = dist(t(eSet2))

rv = rowVars(eSet, na.rm = TRUE)
o = order(rv, decreasing = TRUE)

dists2 = dist(t(eSet[head(o, 6), ]))
# Add based off of top 6

palette(brewer.pal(8, "Dark2"))
hc = hclust(dists2, method = "ward.D2")
dend <- as.dendrogram(hc)
o.dend <- order.dendrogram(dend)
labels(dend) <- colnames(a1)[o.dend]
labels_colors(dend) <- 1*(pData(a1)$summarystage[o.dend] == "superficial") + 2
plot(dend)
