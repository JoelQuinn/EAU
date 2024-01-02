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


data <- readRDS("Documents/EAU_CITE-seq/integrated_data.RDS")

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
# Whole-sample aggregation for PCA ----
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

# Save plot
ggsave("./Data/Plots/SamplePheatmap.png", plot=sampHM, device="png", units="cm",
       dpi=300, width=15, height=10)


# By-cluster aggregation ----
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

splitf

# Create a list where each part of the list is a cluster's 'count matrix',
# then transform so that rows are genes and columns are sample names
agg <- split.data.frame(agg, factor(splitf)) %>%
  lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[[:alnum:]_]+")))

saveRDS(agg, file="Aggregate-for-DESeq2.rds")
sMeta

table(sce$cell_type, sce$orig.ident)

# Cell numbers plot
df <- as.data.frame(table(sce$cell_ids, sce$orig.ident))
df <- rename(df, Cell_Type = Var1, Condition = Var2)


pdf("./Documents/EAU_CITE-seq/Plots/cell-type_numbers.pdf", height=10, width=15)
df %>% ggplot(aes(Cell_Type, Freq)) + 
  geom_col(aes(fill=Condition), position="dodge") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()



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

de_cluster_ids

# Create data frame with sample IDs, cluster names and conditions
gg_df <- data.frame(cluster_id=de_cluster_ids,
                    orig.ident=DEsamples)

gg_df <- left_join(gg_df, sMeta[, c("orig.ident", "condition")])
gg_df

metadata <- gg_df %>% dplyr::select(cluster_id, orig.ident, condition)
metadata

saveRDS(metadata, "pseudobulk_metadata.RDS")

# Run analysis on one cell type - using MG ----
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

#all(rownames(cluster_metadata) == colnames(cluster_counts))
#cluster_counts


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

res
# lfc shrinkage
res <- lfcShrink(dds, coef = 2,
                 res = res)

res

# table of results
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl$log10pvalue <- -log10(res_tbl$pvalue)

res_tbl[res_tbl$gene=="H2-Ab1",]

# Significant genes ----
# set p value cutoff
padj_cutoff <- 0.05

# filter results for significance
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

View(sig_res)

#normalize counts
normalized_counts <- counts(dds, normalized=TRUE)

sig_res
#
sig_res$log10pvalue <- -log10(sig_res$pvalue)

ggplot(data=res_tbl, mapping = aes(x=log2FoldChange, y=log10pvalue)) +
  geom_point() +
  coord_cartesian(xlim=c(-5,10))


# take top 20 significant results
top20 <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)

# plotting top 20 genes (with some wrangling for ggplot2) ----
top20_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var="gene") %>%
  dplyr::filter(gene %in% top20)

gathered_top20 <- top20_norm %>%
  gather(colnames(top20_norm)[2:length(colnames(top20_norm))], key="samplename", value="normalized_counts")

gathered_top20 <- inner_join(sMeta[, c("orig.ident", "condition")], gathered_top20, by=c("orig.ident"="samplename"))

gathered_top20

ggplot(gathered_top20) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = condition), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var="gene") %>%
  dplyr::filter(gene %in% sig_res$gene)


heat_colors <- rev(brewer.pal(10, "RdYlBu"))

rownames(sig_norm) <- sig_norm$gene

cleaned_sig_norm <- sig_norm[!sig_norm$gene %in% rod_genes,]
cleaned_sig_norm["Cd274",]
View(cleaned_sig_norm)

ph <- pheatmap(sig_norm[, 2:length(colnames(sig_norm))], 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("condition", "cluster_id")],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 8, 
         height = 20)




ggsave("./Documents/EAU_CITE-seq/MG_DEseq2_pheatmap.pdf", plot=ph, dpi=300, units="cm", device = "pdf",
       width=10, height=10)

h2 <- cleaned_sig_norm[cleaned_sig_norm$gene %in% grep("H2", rownames(cleaned_sig_norm), value = TRUE),]

h2 <- h2[, c(1, 4, 5, 6, 2, 3)]


phh2 <- pheatmap(h2[, 2:length(colnames(h2))], 
               cluster_rows = T, 
               show_rownames = F,
               annotation = cluster_metadata[, c("condition", "cluster_id")],
               border_color = NA, 
               fontsize = 10, 
               scale = "row", 
               fontsize_row = 8, 
               height = 20, cluster_cols = F) +
  

ggsave("./Data/Plots/MG_DEseq2_H2-pheatmap.png", plot=phh2, dpi=300, units="cm",
       width=15, height=10)


# ----
new_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var="gene") %>%
  pivot_longer(cols=!gene, names_to = "Condition") %>%
  as.data.frame()

new_norm <- rename(new_norm, Sample_ID = Condition)

new_norm$condition <- rep(rep(c("EAU", "Healthy"), c(2, 1)), length(unique(new_norm$gene)))

new_norm <- new_norm[!new_norm$condition == "Resistant", ]


new_norm$condition <- factor(new_norm$condition, levels = c("Healthy", "EAU"))
colnames(new_norm)[4] <- "Condition"

chemokines <- c("Cxcl1","Ccl2", "Cxcl10", "Cxcl12")
mhc <- c("H2-Ab1","H2-Eb1", "H2-Aa")
complement <- c("C3", "C4b", "C1ra")
adh <- c("Vcam1", "Icam1")

pc <- filter(new_norm, gene %in% chemokines) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_bar(aes(fill=Condition), position="dodge", stat="summary", fun="mean", width = 0.7,) +
  geom_point(aes(y=value, group=Condition), position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=22, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=18, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  ylab("Normalized Pseudobulk Counts")

pmhc <- filter(new_norm, gene %in% mhc) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_bar(aes(fill=Condition), position="dodge", stat="summary", fun="mean", width = 0.7,) +
  geom_point(aes(y=value, group=Condition), position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=22, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=18, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

pcomp <- filter(new_norm, gene %in% complement) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_bar(aes(fill=Condition), position="dodge", stat="summary", fun="mean", width = 0.7,) +
  geom_point(aes(y=value, group=Condition), position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=22, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=18, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  coord_cartesian(ylim=c(0, 300))

padh <- filter(new_norm, gene %in% adh) %>%
  ggplot(aes(x=gene, y=value)) +
  geom_bar(aes(fill=Condition), position="dodge", stat="summary", fun="mean", width = 0.7,) +
  geom_point(aes(y=value, group=Condition), position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=22, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=18, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  coord_cartesian(ylim=c(0,300))


pdf("LeadingEdgeMG.pdf", height=5, width=15)
pc + pcomp + padh + pmhc + patchwork::plot_layout(ncol=4)
dev.off()


filter(new_norm, gene %in% "Il6") %>%
  ggplot(aes(x=gene, y=value)) +
  geom_bar(aes(fill=Condition), position="dodge", stat="summary", fun="mean", width = 0.5,) +
  geom_point(aes(y=value, group=Condition), position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=18, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=18, colour="black"),
        legend.text = element_text(size=14)) 

pIgsf <- new_norm[new_norm$gene =="Igsf11",] %>%
  ggplot(aes(x=factor(gene), y=value)) +
  geom_bar(aes(fill=condition), position="dodge", stat="summary", fun="mean", width = 0.5) +
  geom_point(aes(y=value, group=condition), position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=16, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=14, colour="black")) +
  ylab("Normalized Pseudobulk Counts") +
  scale_y_continuous(limits=c(0, 1600), expand=c(0,0)) 

  
  
pmg <-  new_norm[new_norm$gene %in% "Cd274",] %>%
  ggplot(aes(x=factor(gene, levels = goi), y=value)) +
  geom_bar(aes(fill=condition), position="dodge", stat="summary", fun="mean", width = 0.5) +
  geom_point(aes(y=value, group=condition), position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16, colour="black", angle=45, hjust=1, face = "bold.italic"),
        axis.text.y = element_text(size=14, colour="black")) +
  ylab("Normalized Pseudobulk Counts") +
  scale_y_continuous(limits=c(0, 250), expand=c(0,0))

pmg <- pmg + NoLegend()
pmg
fig <- (pmg + theme(axis.title.y = element_blank())) / ((prods | pIgsf) + plot_layout(widths=c(6,1)))


ggsave("./Data/Plots/Rods-Igsf11.pdf", plot=igsf, dpi=600, units="cm", height=10, width=15)

ggsave("./Data/Plots/rods-mg2.pdf", plot=fig, dpi=600, units="cm", height=11, width=22)


colnames(normalized_counts)

write_csv(sig_res, file = "./Data/DESeq2/EAUvsHealthy.csv")

# Genes of interest ----
filter(s_tbl, gene == "H2-Ab1")


#gene sets----
res_tbl <- res_tbl[order(-res_tbl$log2FoldChange),]

gene_list <- res_tbl$log2FoldChange
names(gene_list) <- res_tbl$gene

entrez <- select(org.Mm.eg.db, keys=names(gene_list), columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

names(gene_list) <- entrez$ENTREZID

gene_list

gene_list <- sort(gene_list, decreasing=T)

write.table(gene_list, file="mg_ranked-genes.rnk", sep="\t", quote=F, row.names = T)



load("Documents/EAU_CITE-seq/mouse_H_v5p2.rdata")



pathwaysH <- Mm.H
fgseaRes <- fgsea(pathways = pathwaysH, stats = gene_list, 
                  minSize=15, maxSize=500)


head(fgseaRes[order(padj, -abs(NES)), ], n=20)

fgsea

topPathwaysUP <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDOWN <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUP, rev(topPathwaysDOWN))


fgseaRes[padj < 0.05]

gseaH_plot <- plotGseaTable(pathwaysH[topPathways], stats = gene_list, fgseaRes, gseaParam = 0.5)

pdf("GSEA-H_plot.pdf")
gseaH_plot
dev.off()

fgseaResTidy <- fgseaRes %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

select(org.Mm.eg.db, keys=unlist(fgseaResTidy[8,8]), columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")[,2]


pdf("mg_sig_hallmark_gsea.pdf", height=10, width=10)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj)) +
  scale_color_continuous() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_cowplot()
dev.off()


# Analyse all clusters ----

clusters <- unique(metadata$cluster_id)

clusters
  
for (x in seq(1, length(clusters))) {
  
  # Get cluster name
  cl <- clusters[x]
  
  # Make output directory
  output_dir <- paste("./Documents/EAU_CITE-seq/DESeq2/", cl, sep = "")
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


ggplot(res_tbl, aes(x=log2FoldChange, y=-log(padj))) +
  geom_point()



ncounts <- read.csv(file="./Documents/EAU_CITE-seq/DESeq2/Müller Glia/Müller Glia_norm-counts_EAUvHealthy.csv", header = T)
sigs <- read.csv(file="./Documents/EAU_CITE-seq/DESeq2/Müller Glia/Müller Glia_significant-genes_EAUvHealthy.csv", header = T)





  