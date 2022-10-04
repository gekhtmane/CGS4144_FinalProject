# Make sure BioInfoFile1.csv is in your
# data environment before running this script
library(readr)
geo_data <- read.csv("BioInfoFile1.csv", header = TRUE)
ensembl_ids <- geo_data[, 1]
rownames(geo_data) = ensembl_ids

# convert to Hugo
library("org.Hs.eg.db")
geneSymbols <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
gene_symbols= as.data.frame(as.matrix(geneSymbols))

# store gene types for later use
gene_types <- geo_data[, 9]

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

#t-sne plot using M3C package
install.packages("M3C")
library(M3C)
tsne(geo_data,colvec=c('gold'))


# Part 3
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

#save tsv file
readr::write_tsv(
  deseq_df,
  file.path(
    "DiffExp.tsv" # Replace with a relevant output file name
  )
)

#volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot

# Part 4
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# get statistically significant list of expressed genes from deseq_df
sig_deseq <- deseq_df[deseq_df$pvalue < 0.01, ]

# counts matrix
dds_matrix <- counts(dds[rownames(dds) %in% sig_deseq$Gene])
# get z score values of matrix
z_matrix <- t(apply(dds_matrix, 1, scale))

#FIXME: add better side panel
Heatmap(z_matrix, cluster_rows = T, cluster_columns = T)

# Part 5 - topGO
BiocManager::install("topGO")
library(topGO)

# create named vector with genes and p values
deseq_vector <- setNames(deseq_df$pvalue, deseq_df$Gene)
deseq_vector <- deseq_vector[!is.na(deseq_vector)]

# statistically significant vector
sig_vector <- setNames(sig_deseq$pvalue, sig_deseq$Gene)

# below code is taken from the topGO library documentation
topDiffGenes <- function(pvalue) {
  return(pvalue < 0.01)
}
gene_sel <- topDiffGenes(deseq_vector)

# set annotation to specific db and gene ID
allGO2genes <- annFUN.org(
  whichOnto = "BP",
  feasibleGenes = NULL,
  mapping = "org.Hs.eg.db",
  ID = "ensembl")

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = deseq_vector,
              annot = annFUN.GO2genes,
              GO2genes = allGO2genes,
              geneSel = topDiffGenes)

# now use GOdata to perform enrichment analysis 
# run Fisher test
result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Kolmogorov-Smirnov
result_KS_elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
result_KS_classic <- runTest(GOdata, algorithm = "classic", statistic = "ks")

results_table <- GenTable(GOdata, classicFisher = result_fisher,
                          classicKS = result_KS_classic,
                          elimKS = result_KS_elim, orderBy = "elimKS",
                          ranksOf = "classicFisher", topNodes = 10)

#save csv file
write.table(results_table, file='topGO.csv', quote=FALSE, sep=',')

BiocManager::install("Rgraphviz")
# compare values - code based on topGO R vignette
showSigOfNodes(GOdata, score(result_KS_elim), firstSigNodes = 5, useInfo = 'all')

end()     
 