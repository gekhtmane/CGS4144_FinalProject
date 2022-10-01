#install.packages("htmltools")
#library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
setwd("C:/Users/chidi")

BioInfoFile1_2 <- read.csv('C:/Users/chidi/Downloads/BioInfoFile1.3.csv', header=TRUE)
x <- BioInfoFile1_2[,-(1:9)]

y <- log2(x[,-1]+1) #skip first col which has hugo gene names
z <- apply(X = y, MARGIN = 1, FUN = median)
plot(density(z))

library("DESeq2")
library(ggplot2)

countData <- read.csv('C:/Users/chidi/Downloads/BioInfoFile1.3.csv', header = TRUE, sep = ",")
head(countData)
metaData <- read.csv('C:/Users/chidi/Downloads/SeriesMatrixPatIDSampType.csv', header = TRUE, sep = ",")

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ condition)
dds
##dds <- DESeqDataSetFromMatrix(countData=y, 
##                              colData=metaData, 
##                              design=~dex, tidy = TRUE)
##dds
##dds <- DESeq(dds)
#res <- results(dds)
#head(results(dds, tidy=TRUE))
