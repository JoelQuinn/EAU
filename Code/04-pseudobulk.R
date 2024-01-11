## Pseudobulk analysis with DESeq2

library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(biomaRt)
library(fgsea)
library(msigdbr)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Read in integrated data object
data <- readRDS("path/to/integrated_data.RDS")

# Set up SingleCellExperiment object
counts <- data@assays$RNA@counts
metadata <- data@meta.data

# Create single cell experiment object
sce <- SingleCellExperiment(assays=list(counts=counts),
                            colData=metadata)

# Get groups (original sample and cell type)
groups <- colData(sce)[, c("cell_ids", "orig.ident")]

# Look through object
assays(sce)
dim(counts(sce))

counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))

# Add a metadata column with the overall conditions (i.e. Healthy, Resistant or EAU)
sce$condition <- ifelse(sce$orig.ident == "EAU_1" | sce$orig.ident == "EAU_2", "EAU", sce$orig.ident)


# Generating sample level metadata
# Get number of clusters and their names
clustIDs <- purrr::set_names(levels(sce$cell_ids))
numClust <- length(clustIDs)

# Get number of samples and their names
levels(sce$orig.ident) <- c("Healthy", "EAU_1", "EAU_2")
sampleIDs <- purrr::set_names(levels(sce$orig.ident))
nSamples <- length(sampleIDs)

# Calculate number of cells per sample - NB table reordered to match levels of orig.ident
n_cells <- as.numeric(table(sce$orig.ident)[c(3, 1, 2)])

# Reorder samples to match orig.ident levels
m <- match(sampleIDs, sce$orig.ident)

# Create metadata
sMeta <- data.frame(colData(sce)[m, ], n_cells, row.names=NULL) %>%
  dplyr::select("orig.ident", "condition", "n_cells")


# Count aggregation
# Whole-sample aggregation for PCA
sample_groups <- colData(sce)[, "orig.ident"]

# Aggregate matrix and group by sample
sample_agg <- aggregate.Matrix(t(counts(sce)),
                               groupings=sample_groups,
                               fun="sum") %>% t()

# Create metadata for aggregated matrix
aggMetadata <- sMeta[, c("orig.ident", "condition")]

# Create DESEq object 
samp_dds <- DESeqDataSetFromMatrix(sample_agg,
                                   colData = aggMetadata,
                                   design= ~ condition)

# Normalize counts to log2
samp_rld <- rlog(samp_dds, blind = TRUE)

# Plot PCA by condition
pca_plot <- DESeq2::plotPCA(samp_rld, intgroup="condition")
pca_plot + theme_classic()

# Extract normalized matrix
s_rld_mat <- assay(samp_rld)

# Compute pairwise correlations of each sample
s_rld_cor <- cor(s_rld_mat)

# Make rownames orig.ident for plotting
rownames(aggMetadata) <- aggMetadata$orig.ident

# Heatmap of sample pairwise correlations
sampHM <- pheatmap(s_rld_cor, annotation = aggMetadata[, "condition", drop=F])



# By-cluster aggregation
# subset metadata to only include cluster labels and sample name for aggregation
groups <- colData(sce)[, c("cell_ids", "orig.ident")]

# aggregate across cluster-sample groups
agg <- aggregate.Matrix(t(counts(sce)),
                        groupings=groups, fun="sum")


# output is a matrix of rows = cell type in each sample and columns = genes
agg[1:10, 1:6]

# Create vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(agg),
                          pattern="_",
                          n=2), `[`, 1)

# Create a list where each part of the list is a cluster's 'count matrix',
# then transform so that rows are genes and columns are sample names
agg <- split.data.frame(agg, factor(splitf)) %>%
  lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[[:alnum:]_]+")))

# Get sample names for each cell type cluster

# prepare data.frame for plotting
get_sample_ids <- function(x){
  agg[[x]] %>% colnames()
}

# Get list of sample names repeated for every cluster
DEsamples <- map(1:length(clustIDs), get_sample_ids) %>% unlist()

# get cluster IDs for each of the samples
samples_list <- map(1:length(clustIDs), get_sample_ids)

# For each cluster ID, get the cluster ID and repeat it for each sample
get_cluster_ids <- function(x){
  rep(names(agg)[x],
      each=length(samples_list[[x]]))
}

# Perform function across all cluster IDs
de_cluster_ids <- map(1:length(clustIDs), get_cluster_ids) %>%
  unlist()


# Create data frame with sample IDs, cluster names and conditions
gg_df <- data.frame(cluster_id=de_cluster_ids,
                    orig.ident=DEsamples)

# Add orig.ident and condition 
gg_df <- left_join(gg_df, sMeta[, c("orig.ident", "condition")])

metadata <- gg_df %>% dplyr::select(cluster_id, orig.ident, condition)

# Run analysis on one cell type - using MG (skip to line 280 for function to run on all clusters)

# View available clusters
clusters <- unique(metadata$cluster_id)
clusters
         
# subset MG metadata
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[5]), ]
head(cluster_metadata)

# assign rownames of metadata to be sample id
rownames(cluster_metadata) <- cluster_metadata$orig.ident
head(cluster_metadata)

# subset counts to only MG
counts <- agg[[clusters[5]]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# DESeq2 object creation
dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ condition)

# Transform counts
rld <- rlog(dds, blind = TRUE)

#Plot
p1 <- DESeq2::plotPCA(rld, intgroup="condition")
p1 + theme_classic()

# hierarchical clustering
# Extract rlog matrix and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop=F])

# Run DESeq2 differential expression analysis
dds$condition <- relevel(dds$condition, ref = "Healthy")

dds <- DESeq(dds)

plotDispEsts(dds)

# Output results of Wald test for contrast for EAU vs Healthy
contrast <- c("condition", unique(cluster_metadata$condition)[1], unique(cluster_metadata$condition)[2])
unique(cluster_metadata$condition)

res <- results(dds,
               contrast=contrast,
               alpha = 0.05)

resultsNames(dds)

# lfc shrinkage
res <- lfcShrink(dds, coef = 2,
                 res = res)

# table of results
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl$log10pvalue <- -log10(res_tbl$pvalue)

# Significant genes ----
# set p value cutoff
padj_cutoff <- 0.05

# filter results for significance
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

View(sig_res)
sig_res$log10pvalue <- -log10(sig_res$pvalue)

#normalize counts
normalized_counts <- counts(dds, normalized=TRUE)

# Filter normalized counts by significant results - can be used for plotting
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var="gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Normalized counts dataframe for plotting
new_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var="gene") %>%
  pivot_longer(cols=!gene, names_to = "Condition") %>%
  as.data.frame()

new_norm <- rename(new_norm, Sample_ID = Condition)

new_norm$condition <- rep(rep(c("EAU", "Healthy"), c(2, 1)), length(unique(new_norm$gene)))

new_norm$condition <- factor(new_norm$condition, levels = c("Healthy", "EAU"))
colnames(new_norm)[4] <- "Condition"

# Rank genes for GSEA
res_tbl <- res_tbl[order(-res_tbl$log2FoldChange),]

gene_list <- res_tbl$log2FoldChange
names(gene_list) <- res_tbl$gene

entrez <- select(org.Mm.eg.db, keys=names(gene_list), columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

names(gene_list) <- entrez$ENTREZID

gene_list <- sort(gene_list, decreasing=T)

# Analyse all clusters ----

clusters <- unique(metadata$cluster_id)
  
for (x in seq(1, length(clusters))) {
  
  # Get cluster name
  cl <- clusters[x]
  
  # Make output directory
  output_dir <- paste("path/to/DESeq2_results/", cl, sep = "")
  dir.create(output_dir)
  
  # subset metadata
  cluster_metadata <- metadata[which(metadata$cluster_id == cl), ]
  
  # assign rownames of metadata to be sample id
  rownames(cluster_metadata) <- cluster_metadata$orig.ident
  
  # subset counts to only cluster
  counts <- agg[[cl]]
    
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  # DESeq2 object creation
  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ condition)
  
  # Transform counts
  rld <- rlog(dds, blind = TRUE)
  
  #Plot PCA
  p1 <- DESeq2::plotPCA(rld, intgroup="condition")
  p1 <- p1 + theme_classic()
  
  ggsave(filename=paste(cl, "_PCA.png", sep=""), path=output_dir, plot = p1, 
         device = "png", dpi = 300, width=10, height=15)
  
  
  # hierarchical clustering
  # Extract rlog matrix and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  p2 <- pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop=F])
  
  ggsave(filename=paste(cl, "_CorHeatmap.png", sep=""), path=output_dir, plot = p2, 
         device = "png", dpi = 300, width=15, height=10)
  
  
  # Run DESeq2 differential expression analysis
  dds$condition <- relevel(dds$condition, ref="Healthy")
  
  dds <- DESeq(dds)
  
  # Plot dispersion estimates
  plotDispEsts(dds)
  
  
  ggsave(filename=paste(cl, "_DispEsts.png", sep=""), path=output_dir, 
         device = "png", dpi = 300, width=10, height=15)
  
  
  levels(cluster_metadata$condition) <- c("EAU", "Healthy")
  
  # Set up conditions to compare
  contrast <- c("condition", levels(cluster_metadata$condition)[1], levels(cluster_metadata$condition)[2])
  
  contrast
  
  res <- results(dds,
                 contrast=contrast,
                 alpha = 0.05)
  
  
  resultsNames(dds)
  
  # lfc shrinkage
  res <- lfcShrink(dds,
                   res=res, 
                   coef = 2,
                   type="apeglm")
  
  # table of results
  s_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  
  # set p value cutoff
  padj_cutoff <- 0.05
  
  # filter results for significance
  sig_res <- dplyr::filter(s_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  #normalize counts
  normalized_counts <- counts(dds, normalized=TRUE)
  
  # Rank genes by log2fc for fgsea
  gene_list <- s_tbl$log2FoldChange
  names(gene_list) <- s_tbl$gene
  
  entrez <- select(org.Mm.eg.db, keys=names(gene_list), columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")
  
  names(gene_list) <- entrez$ENTREZID
  
  gene_list
  
  gene_list <- sort(gene_list, decreasing=T)
  
  write.csv(gene_list, file=paste(output_dir, "/", cl, "_ranked-genes_ENTREZ.csv", sep=""))
  write.csv(sig_res, file = paste(output_dir, "/", cl, "_significant-genes_EAUvHealthy.csv", sep=""))
  write.csv(normalized_counts, file=paste(output_dir, "/", cl, "_norm-counts_EAUvHealthy.csv", sep=""))

}




  
