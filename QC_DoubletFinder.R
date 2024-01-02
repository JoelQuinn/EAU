library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(DoubletFinder)
library(DropletUtils)

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


# Run dim reduction and clustering on each sample separately prior to DoubletFinder
samples <- lapply(samples, function(x) SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo")))

samples <- lapply(samples, function(x) RunPCA(x, npcs=75, verbose=F))

ElbowPlot(samples[[1]], ndims=75)

ElbowPlot(samples[[3]], ndims=75)


# Healthy standard processing
healthy <- RunUMAP(samples[[1]], dims = 1:35) %>% 
  FindNeighbors(dims=1:35) %>% 
  FindClusters(resolution=seq(0.3, 1.5, 0.1))

clustree(healthy, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

Idents(healthy) <- "SCT_snn_res.0.8"
hdim <- DimPlot(healthy, reduction="umap", label=TRUE, repel=T) + NoLegend()
hdim

# MG-rod doublet removal
healthy <- FindSubCluster(healthy, cluster=6, graph.name = "SCT_snn", subcluster.name = "MG_doublets")

FindMarkers(healthy, group.by = "MG_doublets", ident.1 = "6_1", ident.2="6_2")

healthy <- subset(healthy, subset = MG_doublets == "6_1", invert=TRUE)
  
# Reprocess
healthy <- SCTransform(healthy, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

healthy <- RunPCA(healthy, npcs = 75)

ElbowPlot(healthy, ndims=75)

healthy <- RunUMAP(healthy, dims = 1:35)


# DoubletFinder
sweep_list <- paramSweep_v3(healthy, PCs=1:35, sct = T)
sweep_stats <- summarizeSweep(sweep_list)
bcmvn <- find.pK(sweep_stats)

bcmvn
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# This section is just trialling different parameter values. Went with nExp of 3% and pK of 0.005.
nExp <- round(0.02*nrow(healthy@meta.data))

healthy <- doubletFinder_v3(healthy, PCs=1:35, pK=0.005, nExp = nExp, sct=T)

healthy[[]]

DimPlot(healthy, group.by = "DF.classifications_0.25_0.005_69")
VlnPlot(healthy, features="nFeature_RNA", group.by = "DF.classifications_0.25_0.005_69")

# Subset data
healthy <- subset(healthy, subset = DF.classifications_0.25_0.005_69 == "Singlet")

# Save as RDS
saveRDS(healthy, file="healthy_SoupX-DF.RDS")




# EAU1
ElbowPlot(samples[[2]], ndims=75)

eau1 <- RunUMAP(samples[[2]], dims = 1:45) %>% 
  FindNeighbors(dims=1:45) %>% 
  FindClusters(resolution=seq(0.3, 1.5, 0.1))

clustree(eau1, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

Idents(eau1) <- "SCT_snn_res.0.8"
DimPlot(eau1, reduction="umap", label=TRUE, repel=T) + NoLegend()


# MG-rod doublet removal
eau1 <- FindSubCluster(eau1, cluster=6, graph.name = "SCT_snn", subcluster.name = "MG_subclusters")
DimPlot(eau1, group.by = "MG_subclusters", label=T) + NoLegend()

FindMarkers(eau1, group.by = "MG_subclusters", ident.1 = "6_1", ident.2=c("6_0", "6_2", "6_3"))

eau1 <- subset(eau1, subset = MG_subclusters == "6_1", invert=TRUE)


# Reprocess
eau1 <- SCTransform(eau1, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

ElbowPlot(eau1, ndims=75)

eau1 <- RunUMAP(eau1, dims = 1:45) %>% 
  FindNeighbors(dims=1:45) %>% 
  FindClusters(resolution=seq(0.3, 1.5, 0.1))

clustree(eau1, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

Idents(eau1) <- "SCT_snn_res.0.8"

DimPlot(eau1, label=T) + NoLegend()

# DoubletFinder
sweep_list <- paramSweep_v3(eau1, PCs=1:45, sct = T)
sweep_stats <- summarizeSweep(sweep_list)
bcmvn <- find.pK(sweep_stats)

bcmvn
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# This section is just trialling different parameter values. Went with nExp of 3% and pK of 0.005.
nExp <- round(0.02*nrow(healthy@meta.data))

eau1 <- doubletFinder_v3(eau1, PCs=1:45, pK=0.02, nExp = nExp, sct=T)

eau1[[]]

DimPlot(eau1, group.by = "DF.classifications_0.25_0.02_68")
VlnPlot(eau1, features="nFeature_RNA", group.by = "DF.classifications_0.25_0.02_68")

# Subset data
eau1 <- subset(eau1, subset = DF.classifications_0.25_0.02_68 == "Singlet")

# Save as RDS
saveRDS(eau1, file="healthy_SoupX-DF.RDS")




# EAU2
ElbowPlot(samples[[3]], ndims=75)

eau2 <- RunUMAP(samples[[3]], dims = 1:50) %>% 
  FindNeighbors(dims=1:50) %>% 
  FindClusters(resolution=seq(0.3, 1.5, 0.1))

clustree(eau2, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

Idents(eau2) <- "SCT_snn_res.1"

DimPlot(eau2, reduction="umap", label=TRUE, repel=T) + NoLegend()


# MG-rod doublet removal
eau2 <- FindSubCluster(eau2, cluster=7, graph.name = "SCT_snn", subcluster.name = "MG_subclusters")
DimPlot(eau2, group.by = "MG_subclusters", label=T) + NoLegend()

FindMarkers(eau2, group.by = "MG_subclusters", ident.1 = "6_2", ident.2=c("6_0", "6_1"))

eau2 <- subset(eau2, subset = MG_subclusters == "7_1", invert=TRUE)


# Reprocess
eau2 <- SCTransform(eau2, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

ElbowPlot(eau2, ndims=75)

eau2 <- RunUMAP(eau2, dims = 1:50) %>% 
  FindNeighbors(dims=1:50) %>% 
  FindClusters(resolution=seq(0.3, 1.5, 0.1))

clustree(eau2, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)

Idents(eau2) <- "SCT_snn_res.0.8"

DimPlot(eau2, label=T) + NoLegend()

# DoubletFinder
sweep_list <- paramSweep_v3(eau2, PCs=1:50, sct = T)
sweep_stats <- summarizeSweep(sweep_list)
bcmvn <- find.pK(sweep_stats)

bcmvn
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# This section is just trialling different parameter values. Went with nExp of 3% and pK of 0.005.
nExp <- round(0.04*nrow(eau2@meta.data))

eau2 <- doubletFinder_v3(eau2, PCs=1:50, pK=0.07, nExp = nExp, sct=T)

eau2[[]]

DimPlot(eau2, group.by = "DF.classifications_0.25_0.07_168")
VlnPlot(eau2, features="nFeature_RNA", group.by = "DF.classifications_0.25_0.07_68")

# Subset data
eau2 <- subset(eau2, subset = DF.classifications_0.25_0.07_68 == "Singlet")

# Save as RDS
saveRDS(eau2, file="EAU2_SoupX-DF.RDS")

