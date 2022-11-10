# reload data to look nicer - alex' code
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


#delete?
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

TopFiveThousandGenes <- geo_data[1:5000,]
TopFiveThousandScale <- scaled_data[1:5000,]


#PAM time!
install.packages(c("cluster", "factoextra"))
library(cluster)
library(factoextra)

FiveThousandPam <- pam(TopFiveThousandGenes, 2)
plot(FiveThousandPam)


#PAM Clustering: 10 genes
TopTenGenes <- geo_data[1:10,]
TopTenScale <- scaled_data[1:10,]
TenPam <- pam(TopTenGenes, 2)
plot(TenPam) #doesn't run, too few genes


#PAM Clustering: 100 genes
TopHundredGenes <- geo_data[1:100,]
TopHundredScale <- scaled_data[1:100,]
HundredPam <- pam(TopHundredGenes, 2)
plot(HundredPam) #doesn't run, too few genes


#PAM Clustering: 1000 genes
TopThousandGenes <- geo_data[1:1000,]
TopThousandScale <- scaled_data[1:1000,]
ThousandPam <- pam(TopThousandGenes, 2)
plot(ThousandPam)


#PAM Clustering: 10000 genes
TopTenThousandGenes <- geo_data[1:10000,]
TopTenThousandScale <- scaled_data[1:10000,]
TenThousandPam <- pam(TopTenThousandGenes, 2)
plot(TenThousandPam)


#PAM alluvial plot
install.packages("ggalluvial")
library("ggalluvial")


#grace alluvial code:
cluster_data <- data.frame(TenPam[[2]]$consensusClass, results100[[2]]$consensusClass,
                           results1000[[2]]$consensusClass, results5000[[2]]$consensusClass, col_data)
colnames(cluster_data) <- c("10Genes", "100Genes", "1000Genes", "5000Genes", "TissueType")
#gg plot
ggplot(as.data.frame(cluster_data),
       aes( y = 128, axis1 = '10Genes', axis2 = '100Genes', axis3 = '1000Genes', axis4 = '5000Genes')) +
  geom_alluvium(aes(fill = TissueType)) +
  geom_stratum(aes(), width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10Genes", "100Genes", "1000Genes", "5000Genes"), expand = c(.05, .05)) +
  ggtitle("ConsensusClusterPlus by number of Genes")



#make a new dataframe with table of all of the PAM clustering results:
#making columns for our data frame
col10primary = table(TenPam$clustering)[[1]]
col10ascites = table(TenPam$clustering)[[2]]
col100primary = table(HundredPam$clustering)[[1]]
col100ascites = table(HundredPam$clustering)[[2]]
col1000primary = table(ThousandPam$clustering)[[1]]
col1000ascites = table(ThousandPam$clustering)[[2]]
col5000primary = table(FiveThousandPam$clustering)[[1]]
col5000ascites = table(FiveThousandPam$clustering)[[2]]
col10000primary = table(TenThousandPam$clustering)[[1]]
col10000ascites = table(TenThousandPam$clustering)[[2]]

levels(col_data$TissueType)

#redo
pamResults <- data.frame(
  c(col10ascites, col10primary),
  c(col100ascites, col100primary),
  c(col1000ascites, col1000primary),
  c(col5000ascites, col5000primary),
  c(col10000ascites, col10000primary),
  levels(col_data$TissueType)
)
pamResults

#make the column names nicer:
names(pamResults)[1] <- "10 Genes"
names(pamResults)[2] <- "100 Genes"
names(pamResults)[3] <- "1000 Genes"
names(pamResults)[4] <- "5000 Genes"
names(pamResults)[5] <- "10000 Genes"
names(pamResults)[6] <- "Tissue Type"

#final table post-manipulation
pamResults


#check if it is properly formatted
is_alluvia_form(as.data.frame(pamResults), axes = 1:6, silent = TRUE)

#hmm KEEP FOR NOW AAA
ggplot(as.data.frame(pamResults), 
       aes(y = 110, axis1 = '10 Genes', axis2 = '100 Genes', axis3 = '1000 Genes', axis4 = '5000 Genes', axis5 = '10000 Genes')) +
  geom_alluvium(aes(fill = `Tissue Type`)) +
  geom_stratum() +
  scale_x_discrete(limits = c("10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"), expand = c(.05, .05)) +
  ggtitle("PAM Clustering by Number of Genes")


#PAM heatmap time wooo
library(ComplexHeatmap)
Heatmap(FiveThousandPam)

library(circlize)
# get statistically significant list of expressed genes from deseq_df
sig_deseq <- deseq_df[deseq_df$pvalue < 0.001, ]

# counts matrix
dds_matrix <- counts(dds[rownames(dds) %in% sig_deseq$Gene])

# get z score values of matrix
z_matrix <- t(apply(dds_matrix, 1, scale))

#Heatmap of matrix, side map
#FiveThousandPam$data
#z_matrix
Heatmap(FiveThousandPam$data, row_km = 5, 
        col = colorRamp2(c(0, 1, 2), c("green", "white", "red")),
        show_column_names = FALSE, 
        show_row_names = FALSE, 
        row_title = NULL, 
        show_row_dend = TRUE)

#col_data$TissueType
Heatmap(FiveThousandPam$data, name = "sample groupings",
        top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
        width = unit(15, "mm"))


#PAM Chi-Square Test
#look at 5000 genes instead of the 110 participants
percAscites = 29/110
percPrimary = 81/110
propAscites = percAscites * 5000
propPrimary = percPrimary * 5000
table(FiveThousandPam$clustering)[[1]]
table(FiveThousandPam$clustering)[[2]]

chisq_df <- data.frame(
  c(propPrimary, propAscites),
  c(table(FiveThousandPam$clustering)[[1]], table(FiveThousandPam$clustering)[[2]])
)
chisq_df

#make the column names nicer:
names(chisq_df)[1] <- "Expected"
names(chisq_df)[2] <- "PAMClustering"

#make the row names nicer:
rownames(chisq_df) <- c("Primary","Ascites")

#adj p-val from chi-square test (2.2e-16)
chisq.test(chisq_df$Expected, chisq_df$PamClustering, correct=FALSE)
p.adjust(2.2e-16, method="BH")


#deleted everything past:
#Graphing Time