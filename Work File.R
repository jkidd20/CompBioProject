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
data(GSE31189_eset)
a1 = GSE31189_eset

## Expressions have negative values. Assuming they are normalized.
rexp = exprs(a1)

