# Make sure BioInfoFile1.csv is in your
# data environment before running this script

# Part 1
library(readr)
geo_data <- read.csv("BioInfoFile1.csv", header = TRUE)

geo_data <- geo_data[,-(1:9)]

#boxplot(log2(x+1))
log_data <- log2(geo_data + 1)
data_medians <- apply(X = log_data, MARGIN = 1, FUN = median)
plot(density(data_medians))

# Part 2 - DESeq2
library(BiocManager)
library(DESeq2)

col_data <- read.csv("SeriesMatrixPatIDSampType.csv", header=TRUE)
# transpose the matrix to match it with the counts data for DESeq2
col_data <- t(col_data)

# make the first row of col_data a header
col_name <- col_data[1, ]
col_data <- data.frame(col_data[2:nrow(col_data), ])
colnames(col_data) <- c(col_name)

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

# transform dataset to create PCA plot
normalize_dds <- vst(dds)
plotPCA(normalize_dds, intgroup = c("TissueType"))

#t-sne plot using M3C package
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("M3C", force = TRUE)
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

deseq_object <- DESeq(dds)
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
    getwd(),
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
end()       
 
