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
library(limma)

eSet = exprs(a1)
tumType = pData(a1)$summarystage
age = pData(a1)$age
gender = pData(a1)$gender
design <- model.matrix(~ tumType)
fit <- lmFit(eSet, design)
fit <- eBayes(fit)
tab <- topTable(fit, n=nrow(eSet), sort.by="none")
