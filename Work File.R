# Playing Around
# To install and load data - Curated Bladder
library(BiocInstaller)
biocLite("curatedBladderData")
library(curatedBladderData)

browseVignettes("curatedBladderData")
data(package="curatedBladderData")

# Link to a good reference paper
# https://www.ncbi.nlm.nih.gov/pubmed/15930339

# Load data
data(GSE1827_eset)
