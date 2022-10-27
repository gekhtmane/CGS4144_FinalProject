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

#PAM PART
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


# video tutorial alluvial plot:
alluvial <- pamResults %>%
  group_by("10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes") %>%
  #summarise(obs = n()) %>%
  ggplot(aes(y = "Tissue Type", axis1 = "10 Genes", axis2 = "100 Genes", axis3 = "1000 Genes", axis4 = "5000 Genes", axis5 = "10000 Genes")) +
  geom_alluvium(aes(fill = "10 Genes"), width = 0,
                knot.pos = 0, reverse = FALSE) +
  geom_stratum(width = 1/12, reverse = FALSE) +
  geom_label(stat = "stratum",
             aes(label = after_stat(stratum), vjust=-3)) +
  scale_x_discrete(limits = c("10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"),
                   expand = c(0.05, 0.05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Alluvial PAM Plot")
alluvial


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

chisq.test(chisq_df$Expected, chisq_df$PamClustering, correct=FALSE)





# Graphing time
library(factoextra)
set.seed(123)
km.res <- kmeans(TopScale, 2, nstart = 25)
fviz_nbclust(TopScale, kmeans, method = "silhouette")
# Visualize kmeans clustering
fviz_cluster(km.res,TopGenes, ellipse.type = "norm", geom = c("point"), main = "5000 Genes")
alluvial_data <- data.frame(km.res$cluster)
# make list for numbers of genes
nums <- list(10, 100, 1000, 10000)
set.seed(123)
for(x in nums){
  TopGenes <- geo_data[1:x,]
  TopScale <- log_data[1:x,]
  TopGenes <- t(TopGenes)
  TopScale <- t(TopScale)
  km.res <- kmeans(TopScale, 2, nstart = 25)

  alluvial_data$x <- km.res$cluster
}
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
