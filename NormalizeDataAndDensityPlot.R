# Make sure BioInfoFile1.csv is in your
# data environment before running this script
library(readr)
geo_data <- read.csv("BioInfoFile1.csv", header = TRUE)

geo_data <- geo_data[,-(1:9)]

#boxplot(log2(x+1))
log_data <- log2(geo_data + 1)
data_medians <- apply(X = log_data, MARGIN = 1, FUN = median)
plot(density(data_medians))

# DESeq2
library(BiocManager)
library(DESeq2)

col_data <- read.csv("SeriesMatrixPatIDSampType.csv", header=TRUE)
# transpose the matrix to match it with the counts data for DESeq2
col_data <- t(col_data)

# make the first row of col_data a header
col_name <- col_data[1, ]
col_data <- data.frame(col_data[2:nrow(col_data), ])
colnames(col_data) <- c(col_name)

all(rownames(col_data) == colnames(geo_data))
geo_data <- geo_data[, rownames(col_data)]
all(rownames(col_data) == colnames(geo_data))

#Formats the data into a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = geo_data, colData = col_data,
                              design = ~ TissueType)
dds <- DESeq(dds)

# transform dataset to create PCA plot
normalize_dds <- vst(dds)
plotPCA(normalize_dds, intgroup = c("TissueType"))

end()     
 