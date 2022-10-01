# First make sure BioInfoFile1 is in your
# data environment before running this script
BioInfoFile1 <- read.csv("C:/Users/surfc/Desktop/BioInformatics/BioInfoFile1.csv", header=TRUE)
x <- BioInfoFile1[,-(1:9)]
#boxplot(log2(x+1))
y <- log2(x+1)
z <- apply(X = y, MARGIN = 1, FUN = median)
plot(density(z))

# DESeq2
library("DESeq2")
ColFilter <- read.csv("C:/Users/surfc/Desktop/BioInformatics/SeriesMatrixPatIDSampType.csv", header=TRUE)
ColFilter <- t(ColFilter)
all(rownames(ColFilter) %in% colnames(x))
x <- x[, rownames(ColFilter)]
all(rownames(ColFilter) == colnames(x)) 
#This Literally just formats the data correctly, using V1 causes and error, Unsure of why this is a problem
dds <- DESeqDataSetFromMatrix(countData = x,
                              colData = ColFilter,
                              design = ~ V1)
dds <- DESeq(dds)
res <- results(dds)
plotCounts(dds,gene=which.min(res$padj), intgroup = "") #DOESNT WORK YET!!!!!!
end()       
 