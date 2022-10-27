
geo_data <- read.csv("BioInfoFile1.csv", header = TRUE)
ensembl_ids <- geo_data[, 1]
rownames(geo_data) = ensembl_ids

# convert to Hugo
library("org.Hs.eg.db")
geneSymbols <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
gene_symbols= as.data.frame(as.matrix(geneSymbols))

geo_data <- geo_data[,-(1:9)]

# DESeq2
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

library(DESeq2)
#Formats the data into a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = geo_data, colData = col_data,
                              design = ~ TissueType)
dds <- DESeq(dds)
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

# NOTE: Deseq2 was already run on dds previously
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
# get statistically significant list of expressed genes from deseq_df
sig_deseq <- deseq_df[deseq_df$pvalue < 0.001, ]

# counts matrix
dds_matrix <- counts(dds[rownames(dds) %in% sig_deseq$Gene])

# start setting up for ConsensusClusterPlus
d <- dds_matrix

# 5000 genes
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:5000],]

# gene median center d, as in tutorial
d = sweep(d,1, apply(d,1,median,na.rm=T))

library(ConsensusClusterPlus)
results5000 <- ConsensusClusterPlus(d, maxK=5, reps=50, pItem=0.8, pFeature=1,
                                 title="5000",clusterAlg="hc",distance="pearson",
                                 seed=1262118388.71279,plot="png")
icl = calcICL(results5000,title="5000",plot="png")

# 10 genes
d <- dds_matrix
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:10],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results10 <- ConsensusClusterPlus(d, maxK=5, reps=50, pItem=0.8, pFeature=1,
                                    title="10",clusterAlg="hc",distance="pearson",
                                    seed=1262118388.71279,plot="png")
icl = calcICL(results10,title="10",plot="png")

# 100
d <- dds_matrix
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:100],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results100 <- ConsensusClusterPlus(d, maxK=5, reps=50, pItem=0.8, pFeature=1,
                                  title="100",clusterAlg="hc",distance="pearson",
                                  seed=1262118388.71279,plot="png")
icl = calcICL(results100,title="100",plot="png")

# 1000
d <- dds_matrix
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:1000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results1000 <- ConsensusClusterPlus(d, maxK=5, reps=50, pItem=0.8, pFeature=1,
                                   title="1000",clusterAlg="hc",distance="pearson",
                                   seed=1262118388.71279,plot="png")
icl = calcICL(results1000,title="1000",plot="png")

# 10000
d <- dds_matrix
mads <- apply(d,1,mad)
d <- d[rev(order(mads))[1:10000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results10000 <- ConsensusClusterPlus(d, maxK=5, reps=50, pItem=0.8, pFeature=1,
                                   title="10000",clusterAlg="hc",distance="pearson",
                                   seed=1262118388.71279,plot="png")
icl = calcICL(results10000,title="10000",plot="png")

#Alluvial Plot
library(ggalluvial)
cluster_data <- data.frame(results10[[2]]$consensusClass, results100[[2]]$consensusClass,
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


# HeatMap
library(ComplexHeatmap)

# chi square
chi_test <- chisq.test(cluster_data$`5000Genes`, cluster_data$TissueType)
p.adjust(chi_test$p.value, method = "BH")
