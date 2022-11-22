# ConsensusClusterPlus
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

library(readr)
geo_data <- read.csv("BioInfoFile1.csv", header = TRUE)
ensembl_ids <- geo_data[, 1]
rownames(geo_data) = ensembl_ids

# convert to Hugo
library("org.Hs.eg.db")
geneSymbols <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
gene_symbols= as.data.frame(as.matrix(geneSymbols))

geo_data <- geo_data[,-(1:9)]

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

# filter to just primary tumor and ascites data
limited <- data.frame(col_data[col_data$TissueType == "Ascites" |
                                 col_data$TissueType == "PrimaryTumour", ], 
                      row.names = rownames(col_data)[col_data$TissueType == "Ascites" |
                                                       col_data$TissueType == "PrimaryTumour"])
colnames(limited) <- c("TissueType")
geo_data <- geo_data[, rownames(limited)]

#Formats the data into a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = geo_data, colData = limited,
                              design = ~ TissueType)
dds <- DESeq(dds)
dds_matrix <- counts(dds[rownames(dds) %in% sig_deseq$Gene])


# now that assignment 2 data is loaded, filter to 5000 most variable genes
d = dds_matrix
mads=apply(d,1,mad)
d = d[rev(order(mads))[1:5000],]

# gene median center d, as indicated in the Cluster tutorial
d = sweep(d,1, apply(d,1,median,na.rm=T))

library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",seed=1262118388.71279, title="5000", plot="png")
icl = calcICL(results,title="5000", plot="png")

# repeat with 10 genes
d = dds_matrix
mads=apply(d,1,mad)
d = d[rev(order(mads))[1:10],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results10 <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                                clusterAlg="hc",distance="pearson",seed=1262118388.71279, title="10", plot="png")
icl10 = calcICL(results,title="10", plot="png")

# repeat with 100 genes
d = dds_matrix
mads=apply(d,1,mad)
d = d[rev(order(mads))[1:100],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results100 <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                                  clusterAlg="hc",distance="pearson",seed=1262118388.71279, title="100", plot="png")
icl10 = calcICL(results,title="100", plot="png")

# repeat with 1000 genes
d = dds_matrix
mads=apply(d,1,mad)
d = d[rev(order(mads))[1:1000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results1000 <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                                  clusterAlg="hc",distance="pearson",seed=1262118388.71279, title="1000", plot="png")
icl1000 = calcICL(results,title="1000", plot="png")

# repeat with 10000 genes
d = dds_matrix
mads=apply(d,1,mad)
d = d[rev(order(mads))[1:10000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results10000 <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                                  clusterAlg="hc",distance="pearson",seed=1262118388.71279, title="10000", plot="png")
icl10000 = calcICL(results,title="10000", plot="png")

# compare across gene levels
cluster_data <- data.frame(results10[[3]]$consensusClass, 
                           results100[[3]]$consensusClass,
                           results1000[[3]]$consensusClass,
                           results[[3]]$consensusClass,
                           limited)
colnames(cluster_data) <- c("10Genes", "100Genes", "1000Genes", "5000Genes", "TissueType");

# attach libraries for alluvial plot
remove.packages("vctrs")
install.packages("vctrs")
install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
install.packages("ggalluvial")
library(tidyselect)
library(ggalluvial)
is_alluvia_form(as.data.frame(cluster_data), axes = 1:5, silent = TRUE)

cluster_data <- to_alluvia_form(as.data.frame(cluster_data), id = TissueType)


ggplot(as.data.frame(cluster_data),
       aes( y = 128, axis1 = '10Genes', axis2 = '100Genes', axis3 = '1000Genes', axis4 = '5000Genes')) +
  geom_alluvium(aes(fill = TissueType)) +
  geom_stratum(aes(), width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10Genes", "100Genes", "1000Genes", "5000Genes"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  ggtitle("ConsensusClusterPlus Clusters by number of Genes")


# chi-square test
test_data <- results[[2]]
test_data <- test_data[[3]]
test_data <- data.frame(test_data, col_data)
chi_test <- chisq.test(test_data$test_data, test_data$TissueType)
p.adjust(chi_test$p.value)
