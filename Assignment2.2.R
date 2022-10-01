#setwd("C:/Users/emgek/Desktop/College/CGS 4144/assignment 2")
x <- read_excel('Assignment2_Edited_Supplementary_File.xlsx')
y <- log2(x[,-1]+1) #skip first col which has hugo gene names
z <- apply(X = y, MARGIN = 1, FUN = median)
plot(density(z))

library("DESeq2")
#plotPCA()

#counts_data
head(x)

#metadata abt tissue type (primary tumor vs ascites)
colData <- read_excel('SeriesMatrixPatIDSampType.xlsx')
colData <- t(colData)

#col names in countsData (x) must match row names in colData (metadata with tissue type)
#skip first col in x, skip first row in colData (must transpose)
all(colnames(x[,-1]) %in% rownames(colData[-1,]))
