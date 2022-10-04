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
                              design = ~ SourceName)
dds <- DESeq(dds)

# transform dataset to create PCA plot
normalize_dds <- vst(dds)

plotPCA(normalize_dds, intgroup = c("Tissue Type"))
#takes multiple gene data, much more than 4 dimensions/variables, and helps
#make a 2d visual of the clustered data. can also reveal which gene is the most
#responsible for clustering the data and tell us how accurate the 2d graph is

install.packages("M3C")
library(M3C)

tsne(geo_data,colvec=c('gold')) 
# we use tsne with t distribution to visualize the 3d representation of the gene
#distribution in a 2d plane

##Tdata <- read.csv("BioInfoFile1.csv", header = TRUE)  data(geo_data)
##umap(geo_data,colvec=c('skyblue'))
##bio.umap <- umap(geo_data)

##bio.umap

##head(bio.umap$layout, 9)
##plot.bio(bio.umap, bio.labels)

install.packages("gprofiler2")
library(gprofiler2)
gostres <- gost(query = c("GO:0005198"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
names(gostres)
head(gostres$result, 3)
##multi_gostres2 <- gost(query = list("chromX" = c("GO:0005198"),
##                                    "chromY" = c("GO:0005829"),
##                                    "chromZ" = c("GO:0070268")), 
##                       multi_query = TRUE)
##names(gostres2)
##head(gostres2$result, 3)
gostplot(gostres, capped = TRUE, interactive = TRUE)
head(gostres$result, 3)
end()       
##names(gostres$meta)
##gostres2 <- gost(query = c("ENSG00000000003", "TSPAN6"), 
##                 organism = "hsapiens", ordered_query = FALSE, 
##                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
##                 measure_underrepresentation = FALSE, evcodes = TRUE, 
##                 user_threshold = 0.05, correction_method = "g_SCS", 
##                 domain_scope = "annotated", custom_bg = NULL, 
##                 numeric_ns = "", sources = NULL)
##head(gostres2$result, 3) 