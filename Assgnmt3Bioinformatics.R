# Get K-means of data

#First we filter 5000 of the most significant genes
library(readr)
library(tidyverse)
library(dplyr)
gene_data <- data.frame(read.csv("BioInfoFile1.csv", row.names = 1, header = TRUE))
col_data <- data.frame(read.csv("SeriesMatrixPatIDSampType.csv", header=TRUE))
geo_data <- gene_data[,-(1:8)]
scaled_data <- scale(geo_data)

# transpose the matrix to match it with the counts data for DESeq2
col_data <- t(col_data)
# make the first row of col_data a header
col_name <- col_data[1, ]
col_data <- data.frame(col_data[2:nrow(col_data), ])
colnames(col_data) <- c(col_name)

col_data <- col_data %>% filter(grepl('PrimaryTumour|Ascites',TissueType))
geo_data <- t(geo_data)
tokeep <- c(rownames(col_data))
geo_data <- geo_data[rownames(geo_data) %in% tokeep,]
geo_data <- t(geo_data)



alluvial_data <- read.csv("BioInfoFile1.csv", header = TRUE)
alluvial_data <- alluvial_data[,-(2:9)]

log_data <- log2(geo_data + 1)

#DIFFSEQ
library(BiocManager)
library(DESeq2)


#check if in order
all(rownames(col_data) == colnames(geo_data))
#if false:
geo_data <- geo_data[, rownames(col_data)]
all(rownames(col_data) == colnames(geo_data)) #should now return true

col_data$TissueType <- factor(col_data$TissueType)
#Formats the data into a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = geo_data, colData = col_data,
                              design =~ TissueType )
dds <- DESeq(dds)
normalize_dds <- vst(dds)

# these next lines are from the Alex' Lemonade tutorial
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)

#set seed
set.seed(12345)

#metadata file is labelled col_data
#actual data is labelled geo_data (expression_df)

deseq_object <- dds
deseq_results <- results(deseq_object)

#possible replacement for lfcShrink() line in tutorial:
res_lfc <- lfcShrink(deseq_object, coef=2, type="apeglm")
head(res_lfc) #prints table

# filtering/cleaning up using TidyVerse
deseq_df <- res_lfc %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
head(deseq_df)

# Now we filter

# make the first row of col_data a header
geo_data <- geo_data[deseq_df$Gene,]
scaled_data <- scaled_data[deseq_df$Gene,]

TopGenes <- geo_data[1:5000,]
TopScale <- scaled_data[1:5000,]
# Graphing time
library(factoextra)
set.seed(123)
km.res <- kmeans(TopScale, 10, nstart = 25)

# Visualize kmeans clustering
fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main = "5000 Genes")
#After that we run our k-means clustering and graph it
sort(table(km.res$cluster))
km.res$cluster

#Alluvial Plot
library(ggalluvial)
head(as.data.frame(TopGenes))
head(as.data.frame(col_data))

TopGenes$GeneID <- row.names(TopGenes)
TopGenes$Cluster <- km.res$cluster

ggplot(as.data.frame(as.data.frame(TopGenes)),
       aes( axis1 = GeneID, axis2 = Cluster)) +
  geom_alluvium(aes(fill = Admit), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Gene Cluster Membership")
