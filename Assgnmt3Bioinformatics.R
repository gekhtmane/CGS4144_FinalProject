# Get K-means of data

#First we filter 5000 of the most significant genes
library(readr)
library(tidyverse)
library(dplyr)
gene_data <- data.frame(read.csv("BioInfoFile1.csv", row.names = 1, header = TRUE))
col_data <- data.frame(read.csv("SeriesMatrixPatIDSampType.csv", header=TRUE))
geo_data <- gene_data[,-(1:8)]

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

# Order data
geo_data <- geo_data[deseq_df$Gene,]
log_data <- log_data[deseq_df$Gene,]

# For 5000 genes
TopGenes <- geo_data[1:5000,]
TopScale <- log_data[1:5000,]
TopGenes <- t(TopGenes)
TopScale <- t(TopScale)

Top5000Genes <- geo_data[1:5000,]
Top5000Scale <- log_data[1:5000,]
Top5000Genes <- t(Top5000Genes)
Top5000Scale <- t(Top5000Scale)

# Graphing time
library(factoextra)
set.seed(123)
km.res <- kmeans(TopScale, 2, nstart = 25)
fviz_nbclust(TopScale, kmeans, method = "silhouette")
# Visualize kmeans clustering
fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main = "5000 Genes")

alluvial_data <- data.frame(km.res$cluster)
alluvial_data <- alluvial_data %>% rename('km.res.cluster' = '5000Genes')

TissueType <- c("PrimaryTumour", "Ascites")
Freq <- c(81,29)
Test <- c("ActualCount", "ActualCount")
a_data <- data.frame(TissueType, Test, Freq)

TissueType <- c("PrimaryTumour", "Ascites")
Freq <- c(km.res$size[1],km.res$size[2])
Test <- c("5000Genes", "5000Genes")
a_data <- rbind(a_data, data.frame(TissueType, Test, Freq))
# make list for numbers of genes
# 10 genes
  set.seed(123)
  TopGenes <- geo_data[1:10,]
  TopScale <- log_data[1:10,]
  TopGenes <- t(TopGenes)
  TopScale <- t(TopScale)
  km.res <- kmeans(TopScale, 2, nstart = 25)
  new <- km.res$cluster
  new <- data.frame(new)
  alluvial_data <- cbind(alluvial_data, new)
  alluvial_data <- alluvial_data %>% rename(new = '10Genes')
  TissueType <- c("PrimaryTumour", "Ascites")
  Freq <- c(km.res$size[1],km.res$size[2])
  Test <- c("10Genes", "10Genes")
  a_data <- rbind(a_data, data.frame(TissueType, Test, Freq))
  fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main ='10Genes')
  
# 100 genes
  set.seed(123)
  TopGenes <- geo_data[1:100,]
  TopScale <- log_data[1:100,]
  TopGenes <- t(TopGenes)
  TopScale <- t(TopScale)
  km.res <- kmeans(TopScale, 2, nstart = 25)
  new <- km.res$cluster
  new <- data.frame(new)
  alluvial_data <- cbind(alluvial_data, new)
  alluvial_data <- alluvial_data %>% rename(new = '100Genes')
  TissueType <- c("PrimaryTumour", "Ascites")
  Freq <- c(km.res$size[1],km.res$size[2])
  Test <- c("100Genes", "100Genes")
  a_data <- rbind(a_data, data.frame(TissueType, Test, Freq))
  fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main ='100Genes')
  
# 1000 genes
  set.seed(123)
  TopGenes <- geo_data[1:1000,]
  TopScale <- log_data[1:1000,]
  TopGenes <- t(TopGenes)
  TopScale <- t(TopScale)
  km.res <- kmeans(TopScale, 2, nstart = 25)
  new <- km.res$cluster
  new <- data.frame(new)
  alluvial_data <- cbind(alluvial_data, new)
  alluvial_data <- alluvial_data %>% rename(new = '1000Genes')
  TissueType <- c("PrimaryTumour", "Ascites")
  Freq <- c(km.res$size[1],km.res$size[2])
  Test <- c("1000Genes", "1000Genes")
  a_data <- rbind(a_data, data.frame(TissueType, Test, Freq))
  fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main ='1000Genes')
  
# 10000 genes
  set.seed(123)
  TopGenes <- geo_data[1:10000,]
  TopScale <- log_data[1:10000,]
  TopGenes <- t(TopGenes)
  TopScale <- t(TopScale)
  km.res <- kmeans(TopScale, 2, nstart = 25)
  new <- km.res$cluster
  new <- data.frame(new)
  alluvial_data <- cbind(alluvial_data, new)
  alluvial_data <- alluvial_data %>% rename(new = '10000Genes')
  TissueType <- c("PrimaryTumour", "Ascites")
  Freq <- c(km.res$size[1],km.res$size[2])
  Test <- c("10000Genes", "10000Genes")
  a_data <- rbind(a_data, data.frame(TissueType, Test, Freq))
  fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main ='10000Genes')
  

#Alluvial Plot
library(ggalluvial)

new <- data.frame(col_data$TissueType)
alluvial_data <- cbind(alluvial_data, new)
alluvial_data <- alluvial_data %>% rename("col_data.TissueType" = "TissueType")

alluvial_data$Sample <- row.names(alluvial_data)

is_alluvia_form(as.data.frame(alluvial_data))
is_alluvia_form(as.data.frame(a_data))
head(as.data.frame(alluvial_data))

#gg plot
ggplot(data = a_data,
       aes(y = Freq, axis1 = 'Test', axis2 = "TissueType"))+
  geom_alluvium(aes(fill = TissueType),width = 1/12) +
  geom_stratum(width = 1/4, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("TissueType"), expand = c(.05, .05)) +
  xlab("Group Membership")
  scale_fill_brewer(type = "qual", palette = "Set2") +
  ggtitle("Gene Cluster Membership")



# Heatmap 5000 genes
library(ComplexHeatmap)
library(circlize)
# get statistically significant list of expressed genes from deseq_df
sig_deseq <- deseq_df[deseq_df$pvalue < 0.001, ]

# counts matrix
dds_matrix <- counts(dds[rownames(dds) %in% sig_deseq$Gene])

# get z score values of matrix
z_matrix <- t(apply(Top5000Genes, 1, scale))

#Heatmap of matrix, side map
Heatmap(z_matrix, row_km = 5, 
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        show_column_names = FALSE,show_row_names = FALSE, row_title = NULL, show_row_dend = TRUE)

Heatmap(col_data$TissueType, name = "sample groupings",
        top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
        width = unit(15, "mm"))

#Chi-Squared Test
chisq.test(alluvial_data$`5000Genes`, col_data$TissueType, correct=FALSE)
chisq.test(alluvial_data$`10Genes`, alluvial_data$`100Genes`, correct=FALSE)
chisq.test(alluvial_data$`10Genes`, alluvial_data$`1000Genes`, correct=FALSE)
chisq.test(alluvial_data$`10Genes`, alluvial_data$`10000Genes`, correct=FALSE)
chisq.test(alluvial_data$`100Genes`, alluvial_data$`1000Genes`, correct=FALSE)
chisq.test(alluvial_data$`100Genes`, alluvial_data$`10000Genes`, correct=FALSE)
chisq.test(alluvial_data$`1000Genes`, alluvial_data$`10000Genes`, correct=FALSE)
chi <- c(102.27, 88.315, 55.495, 94.988, 59.688, 69.122)
pval <- c(2.2e-16, 2.2e-16, 9.367e-14, 2.2e-16, 1.111e-14, 2.2e-16)
p.adjust(pval,method = "BH")
padj <- c(3.3000e-16, 3.3000e-16, 9.3670e-14, 3.3000e-16, 1.3332e-14, 3.3000e-16)
chisqrtbl <- data.frame(chi)
chisqrtbl <- cbind(pval)
chisqrtbl <- data.frame(chi)

