library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(DoubletFinder)

setwd("./Documents/EAU_CITE-seq")

# read files ----
hgx <- Read10X("./SoupX_TY73")
res1gx <- Read10X("SoupX_TY74")
res2gx <- Read10X("SoupX_TY75")
eau1gx <- Read10X("SoupX_TY76")
eau2gx <- Read10X("SoupX_TY77")

# Healthy ----
h <- CreateSeuratObject(hgx, project="Healthy", min.cells = 3, min.features = 200)
# Resistant 1 ----
r1 <- CreateSeuratObject(res1gx, project="Resistant_1", min.cells = 3, min.features = 200)
# Res 2 ----
r2 <- CreateSeuratObject(res2gx, project="Resistant_2", min.cells = 3, min.features = 200)
# EAU 1 ----
e1 <- CreateSeuratObject(eau1gx, project="EAU_1", min.cells = 3, min.features = 200)
# EAU2 ----
e2 <- CreateSeuratObject(eau2gx, project="EAU_2", min.cells = 3, min.features = 200)


# QC metric function
qc_calc <- function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern="mt-")
  x[["percent.ribo"]] <- PercentageFeatureSet(x, pattern="Rp[sl]")
  return(x)
}

h <- qc_calc(h)
r1 <- qc_calc(r1)
r2 <- qc_calc(r2)
e1 <- qc_calc(e1)
e2 <- qc_calc(e2)

samples <- c(h, r1, r2, e1, e2)

# QC plots
qc_scatters <- lapply(samples, function(x){
  p1 <- FeatureScatter(x, feature1="nCount_RNA", feature2="nFeature_RNA", pt.size = 0.01) + NoLegend() +
    scale_x_log10() + theme(axis.text.x = element_text(angle = 45,
                                                       hjust=1))
  p2 <- FeatureScatter(x, feature1="nCount_RNA", feature2="percent.mt", pt.size = 0.01) + NoLegend() +
    scale_x_log10() + theme(axis.text.x = element_text(angle = 45,
                                                       hjust=1))
  p3 <- FeatureScatter(x, feature1="nCount_RNA", feature2="percent.ribo", pt.size = 0.01) + 
    scale_x_log10() + theme(axis.text.x = element_text(angle = 45,
                                                       hjust=1))
  wrap_plots(p1, p2, p3)
})

ggsave(filename = "Healthy_QC_scatter.png", path="QC_plots", plot = qc_scatters[[1]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "Resistant-1_QC_scatter.png", path="QC_plots", plot = qc_scatters[[2]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "Resistant-2_QC_scatter.png", path="QC_plots", plot = qc_scatters[[3]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "EAU-1_QC_scatter.png", path="QC_plots", plot = qc_scatters[[4]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "EAU-2_QC_scatter.png", path="QC_plots", plot = qc_scatters[[5]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)



qcMets <- c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.ribo")

qc_violins <- lapply(samples, function(x) {
  v1 <- VlnPlot(x, features = qcMets[1], pt.size = 0.01) +
    scale_y_log10() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    NoLegend()
  
  v2 <- VlnPlot(x, features = qcMets[2], pt.size = 0.01) +
    scale_y_log10() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    NoLegend()
  
  v3 <- VlnPlot(x, features = qcMets[3], pt.size = 0.01) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    NoLegend()
  
  v4 <- VlnPlot(x, features = qcMets[4], pt.size = 0.01) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  wrap_plots(v1, v2, v3, v4, ncol=4)
})

qc_violins[[4]]

ggsave(filename = "Healthy_QC_violins.png", path="QC_plots", plot = qc_violins[[1]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "Resistant-1_QC_violins.png", path="QC_plots", plot = qc_violins[[2]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "Resistant-2_QC_violins.png", path="QC_plots", plot = qc_violins[[3]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "EAU-1_QC_violins.png", path="QC_plots", plot = qc_violins[[4]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)

ggsave(filename = "EAU-2_QC_violins.png", path="QC_plots", plot = qc_violins[[5]], device = "png", units = "cm", 
       width = 20, height=10, dpi=300)


# Subset based on QC metrics
h <- subset(h, 
            subset = nFeature_RNA > 400 & nFeature_RNA < 4500 & nCount_RNA < 10000 & nCount_RNA > 800 & percent.mt < 15)

r1 <- subset(r1, 
             subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA < 15000 & nCount_RNA > 800 & percent.mt < 15)
r2 <- subset(r2, 
             subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA < 15000 & nCount_RNA > 800 & percent.mt < 15)
e1 <- subset(e1, 
             subset = nFeature_RNA > 400 & nFeature_RNA < 4500 & nCount_RNA < 10000 & nCount_RNA > 800 & percent.mt < 15)
e2 <- subset(e2, 
             subset = nFeature_RNA > 400 & nFeature_RNA < 4500 & nCount_RNA < 10000 & nCount_RNA > 800 & percent.mt < 15)

# Make list of samples to pass through functions
samples <- c(h, e1, e2)



# Perform SCTransform and integration ----
samples <- lapply(samples, function(x) SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo")))
features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 3000)
samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list=samples, normalization.method="SCT", anchor.features = features)
combined_int <- IntegrateData(anchorset=anchors, normalization.method="SCT")

# Dimensionality reduction and clustering ----
combined_int <- RunPCA(combined_int, npcs = 75)
elb <- ElbowPlot(combined_int, ndims=75)
elb

ggsave("Elbow-Plot_default-SCTrandform.png", path="QC_plots", plot=elb, device = "png", units = "cm", 
       width = 15, height=10, dpi=300)

combined_int <- RunUMAP(combined_int, dims = 1:50)
combined_int <- FindNeighbors(combined_int, dims = 1:50)

# Use several resolutions for input into clustree
combined_int <- FindClusters(combined_int, resolution=seq(0.3, 1.5, 0.1))

# Visualize with clustree to determine reasonable cluster resolution
clustree <- clustree(combined_int, prefix="integrated_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

clustree
ggsave("clustree_default.png", path="Integrated_Plots", plot=clustree, device="png", units="cm",
       width=15, height=10, dpi=300)

Idents(combined_int) <- "integrated_snn_res.0.8"


# UMAP plot with unlabelled clusters
dplot1g <- DimPlot(combined_int, reduction="umap", label=TRUE, repel=T) + NoLegend()
dplot1g


saveRDS(combined_int, file = "combined_int.rds")

# Marker genes
DefaultAssay(combined_int) <- "RNA"

combined_int <- NormalizeData(combined_int) %>%
  FindVariableFeatures()

genes <- rownames(combined_int)
combined_int <- ScaleData(combined_int, features = genes)


markers <- FindAllMarkers(combined_int, only.pos=TRUE, min.pct = 0.25)


top10 <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top5 <- markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

write.csv(top10, "Top10_markers.csv")

subset(combined_int, ident=5)



ggsave("filtered-min-ribo2_regressed-mt-ribo_UMAP_40pcs.png", plot=dplot1g, path="Integrated_Plots", device= "png", units = "cm", 
       width = 15, height=10, dpi=300)


# UMAP plots coloured by sample for checking integration
dplot2 <- DimPlot(combined_int, label=F, group.by = "orig.ident", pt.size=0.01)
dplot2

ggsave("default_UMAP_by-ident.png", plot=dplot2, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


dplot3 <- DimPlot(combined_int, label=F, group.by = "orig.ident", split.by="orig.ident", pt.size=0.01)
dplot3

ggsave("default_UMAP_split-ident.png", plot=dplot2, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


resample <- SplitObject(combined_int, split.by="orig.ident")

resample <- lapply(resample, function(x) SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo")))
features <- SelectIntegrationFeatures(object.list = resample, nfeatures = 3000)
resample <- PrepSCTIntegration(object.list = resample, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list=resample, normalization.method="SCT", anchor.features = features)
combined_int <- IntegrateData(anchorset=anchors, normalization.method="SCT")


resample[2]

combined_int <- subset(combined_int, subset = percent.ribo > 2)



# Plotting QC metrics on UMAP
count_umap <- FeaturePlot(combined_int, features="nCount_RNA", label = T, repel=T)

ggsave("count_umap.png", plot=count_umap, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


feat_umap <- FeaturePlot(combined_int, features="nFeature_RNA", label = T, repel=T)

ggsave("feat_umap.png", plot=feat_umap, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


mt_umap <- FeaturePlot(combined_int, features="percent.mt", label = T, repel=T)

ggsave("mt_umap.png", plot=mt_umap, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


ribo_umap <- FeaturePlot(combined_int, features="percent.ribo", label = T, repel=T, max.cutoff = 20)
ribo_umap
ggsave("ribo_umap.png", plot=ribo_umap, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


# pca plots
pca_clusters <- DimPlot(combined_int, reduction = "pca", label=T, repel=T) + NoLegend()

ggsave("clusters_pca.png", plot=pca_clusters, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


pca_idents <- DimPlot(combined_int, reduction = "pca", group.by="orig.ident", split.by="orig.ident")

ggsave("idents_pca.png", plot=pca_idents, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


pca_mt <- FeaturePlot(combined_int, reduction="pca", features="percent.mt", max.cutoff = 20, order = TRUE, label=T, repel=T)

ggsave("mt_pca.png", plot=pca_mt, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)


pca_ribo <- FeaturePlot(combined_int, reduction="pca", features="percent.ribo", max.cutoff = 20, order = TRUE, label=T, repel=T)

ggsave("ribo_pca.png", plot=pca_ribo, path="Integrated_Plots", device= "png", units = "cm", 
       width = 20, height=10, dpi=300)




samples <- lapply(samples, function(x){
  x  <- NormalizeData(x) %>%
    FindVariableFeatures()
})

feats <- SelectIntegrationFeatures(samples)
anchs <- FindIntegrationAnchors(samples, anchor.features = feats)
norm_int <- IntegrateData(anchs)

RunPCA(norm_int, npcs = 75)


DefaultAssay(norm_int) <- "integrated"

norm_int <- ScaleData(norm_int) %>%
  RunPCA()

ElbowPlot(norm_int, ndims = 50)

norm_int <- RunUMAP(norm_int, dims = 1:30)
norm_int <- FindNeighbors(norm_int, reduction = "pca", dims = 1:30)

norm_int <- FindClusters(norm_int, resolution=seq(0.3, 1.5, 0.1))

# Visualize with clustree to determine reasonable cluster resolution
clustree <- clustree(norm_int, prefix="integrated_snn_res.",
                     node_colour="sc3_stability", 
                     node_size_range=c(4,10), 
                     edge_arrow=FALSE, 
                     show_axis=TRUE)


clustree

Idents(norm_int) <- "integrated_snn_res.0.8"

DimPlot(norm_int, label=T, repel=T) + NoLegend()






crud <- subset(combined_int, idents=c(1, 15, 14, 16))
crud

ggplot(data = crud[[]]) + geom_bar(aes(fill=orig.ident, x=integrated_snn_res.0.8), position="dodge")




