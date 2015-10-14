# Analysis of Lawson's cancer stem cell data
# For details on getting data with R go to https://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf
# For details on GEO acronyms go to http://www.ncbi.nlm.nih.gov/geo/info/overview.html

setwd("~/Documents/Graduate School/Research/single_cell_dataset_analysis/Lawson/")
library(GEOquery)
library(limma)
library(Biobase)

# -----------------------------------------------------------------------------
# Make sure GSEMatrix = FALSE
gse <- getGEO("GSE70555",GSEMatrix = FALSE)

# Take a peak at some initial information
head(Meta(gse))
names(GSMList(gse)) # samples in series
names(GPLList(gse)) # platform for series

# Lets look at the first sample
GSMList(gse)[[1]]
# -----------------------------------------------------------------------------

# Load series into matrix
gse <- getGEO("GSE70555",GSEMatrix = TRUE)
if (length(gse) > 1) idx <- grep("GPL20665", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]
X <- exprs(gse)
X[is.na(X)] <- 0

# Save to files
write.table(X,'./expressions.txt',sep="\t",row.names = FALSE,col.names = FALSE)  # Expression matrix
write(colnames(X), file='./cells.txt') # Cell names
write(rownames(X), file='./genes.txt') # Gene names
