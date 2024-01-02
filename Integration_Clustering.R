library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(DoubletFinder)
library(DropletUtils)
library(biomaRt)
library(shadowtext)

healthy <- readRDS("healthy_SoupX-DF.RDS")
eau1 <- readRDS("EAU1_SoupX-DF.RDS")
eau2 <- readRDS("EAU2_SoupX-DF.RDS")




# Make list of samples to pass through functions
samples <- c(healthy, eau1, eau2)

# Perform SCTransform and integration ----
features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 3000)
samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list=samples, normalization.method="SCT", anchor.features = features)
data <- IntegrateData(anchorset=anchors, normalization.method="SCT")

# Dimensionality reduction and clustering ----
data <- RunPCA(data, npcs = 75)
elb <- ElbowPlot(data, ndims=75)
elb

data <- RunUMAP(data, dims = 1:30) %>% FindNeighbors(dims = 1:30)

# Use several resolutions for input into clustree
data <- FindClusters(data, resolution=seq(0.3, 1.5, 0.1))

# Visualize with clustree to determine reasonable cluster resolution
clustree <- clustree(data, prefix="integrated_snn_res.",
                     node_colour="sc3_stability", 
                     node_size_range=c(4,10), 
                     edge_arrow=FALSE, 
                     show_axis=TRUE)

clustree

# Set cluster resolution
Idents(data) <- "integrated_snn_res.0.6"

# Visualise clustering
p1 <- DimPlot(data, label=T) + NoLegend()
p1
# Visualise integration
DimPlot(data, label=T, group.by = "orig.ident")


# Cluster annotation
ids <- c("Rods",
         "Rods",
         "Rods",
         "Rods",
         "Rod Bipolar Cells",
         "Cones",
         "Müller Glia",
         "Immune Cells",
         "Cones",
         "Cone Bipolar Cells",
         "Immune Cells",
         "Müller Glia",
         "Cones",
         "Cone Bipolar Cells",
         "Amacrine Cells",
         "Cone Bipolar Cells",
         "Cones",
         "Rod Bipolar Cells",
         "Cone Bipolar Cells",
         "Amacrine Cells",
         "Amacrine Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Junk",
         "RPE")

names(ids) <- levels(data)
data <- RenameIdents(data, ids)

data[["cell_ids"]] <- Idents(data)

data <- subset(data, idents="Junk", invert=TRUE)



# Re-split samples to pass through functions
samples <- SplitObject(data, split.by = "orig.ident")

# Perform SCTransform and integration ----
samples <- lapply(samples, function(x) SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo")))

features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 3000)
samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list=samples, normalization.method="SCT", anchor.features = features)
data <- IntegrateData(anchorset=anchors, normalization.method="SCT")

# Dim reduction and clustering
data <- RunPCA(data, npcs = 75)
elb <- ElbowPlot(data, ndims=75)
elb

data <- RunUMAP(data, dims = 1:40) %>% FindNeighbors(dims = 1:40)

# Use several resolutions for input into clustree
data <- FindClusters(data, resolution=seq(0.3, 1.5, 0.1))

# Visualize with clustree to determine reasonable cluster resolution
clustree <- clustree(data, prefix="integrated_snn_res.",
                     node_colour="sc3_stability", 
                     node_size_range=c(4,10), 
                     edge_arrow=FALSE, 
                     show_axis=TRUE)

clustree

# Set cluster resolution
Idents(data) <- "integrated_snn_res.0.7"

# Visualise clustering
p1 <- DimPlot(data, label=T) + NoLegend()
p1
# Visualise integration
DimPlot(data, group.by = "orig.ident", split.by = "orig.ident")

DefaultAssay(data) <- "RNA"
data <- NormalizeData(data) %>% ScaleData()


ids <- c("Rods",
         "Rods",
         "Rods",
         "Rods",
         "Rod Bipolar Cells",
         "Cones",
         "Müller Glia",
         "Immune Cells",
         "Cones",
         "Cone Bipolar Cells",
         "Cones",
         "Müller Glia",
         "Immune Cells",
         "Cones",
         "Amacrine Cells",
         "Cone Bipolar Cells",
         "Rod Bipolar Cells",
         "Cone Bipolar Cells",
         "Amacrine Cells",
         "Amacrine Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Cone Bipolar Cells",
         "Rods",
         "RPE",
         "Immune Cells")

names(ids) <- levels(data)
data <- RenameIdents(data, ids)

data[["cell_ids"]] <- Idents(data)

DimPlot(data, label=T) + NoLegend()

pdf(file="./Documents/EAU_CITE-seq/UMAP-by-sample.pdf", width = 8, height = 4)
DimPlot(data, split.by = "orig.ident")
dev.off()


# Immune Cells ----
ics <- subset(data, idents="Immune Cells")

DefaultAssay(ics) <- "RNA"
ics <- SCTransform(ics, method="glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

ics <- RunPCA(ics, npcs = 50)
ElbowPlot(ics, ndims = 50)

DefaultAssay(ics) <- "SCT"
ics <- RunUMAP(ics, dims = 1:39)
ics <- FindNeighbors(ics, dims = 1:39)
ics <- FindClusters(ics, resolution = 1.3)

#clustree(ics, prefix="SCT_snn_res.",
 #        node_colour="sc3_stability", 
  #       node_size_range=c(4,10), 
   #      edge_arrow=FALSE, 
    #     show_axis=TRUE)


Idents(ics) <- "SCT_snn_res.1.3"

DimPlot(ics, label=T) + NoLegend()

DefaultAssay(ics) <- "RNA"
ics <- NormalizeData(ics) %>% ScaleData()


ic_markers <- FindAllMarkers(ics, min.pct = 0.25, only.pos = T) %>% 
  group_by(cluster) %>% slice_max(n=10, order_by = avg_log2FC)


View(ic_markers)

ic_ids <- c("Microglia",
            "Naive CD4+ T cells",
            "Rogdi+ MPs",
            "CD8+ T cells",
            "Monocytes",
            "Tregs",
            "Th1 Cells",
            "NK cells",
            "Neutrophils",
            "pDCs",
            "Th17/gd T cells")

names(ic_ids) <- levels(ics)
ics <- RenameIdents(ics, ic_ids)

DimPlot(ics, label=T) + NoLegend()

table(Idents(ics))

#Removing Rogdi MPs and replotting
ics <- subset(ics, idents = "Rogdi+ MPs", invert=T)

DefaultAssay(ics) <- "RNA"
ics <- SCTransform(ics, method="glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

ics <- RunPCA(ics, npcs = 50)
ElbowPlot(ics, ndims = 50)

DefaultAssay(ics) <- "SCT"
ics <- RunUMAP(ics, dims = 1:20)
ics <- FindNeighbors(ics, dims = 1:20)
ics <- FindClusters(ics, resolution = seq(0.3, 1.5,0.1))

clustree(ics, prefix="SCT_snn_res.",
        node_colour="sc3_stability",
        node_size_range=c(4,10),
        edge_arrow=FALSE,
        show_axis=TRUE)

Idents(ics) <- "SCT_snn_res.1.2"

DimPlot(ics, label=T) + NoLegend()

DefaultAssay(ics) <- "RNA"
ics <- NormalizeData(ics) %>% ScaleData()


ic_markers <- FindAllMarkers(ics, min.pct = 0.25, only.pos = T) %>% 
  group_by(cluster) %>% slice_max(n=10, order_by = avg_log2FC)

View(ic_markers)

ic_ids <- c("Microglia",
            "Naive CD4+ T cells",
            "CD8+ T cells",
            "Monocytes",
            "Tregs",
            "Th1 Cells",
            "NK cells",
            "Neutrophils",
            "Th17/gd T cells",
            "pDCs")

names(ic_ids) <- levels(ics)
ics <- RenameIdents(ics, ic_ids)

pdf(file="./Documents/EAU_CITE-seq/ic-UMAP.pdf")
DimPlot(ics, label=T) + NoLegend()
dev.off()

table(Idents(ics))

#Plotting marker genes
feats <- c("Hexb", "P2ry12",
           "Cd3d", "Cd4", "S1pr1",
           "Cd8a", "Gzmk",
           "Fn1", "Cd14",
           "Foxp3",
           "Ifng", "Tbx21",
           "Klrb1c",
           "Ncr1",
           "Csf3r", "S100a8",
           "Il17a", "Rorc", "Tcrg-C1",
           "Klk1", "Tcf4")


pdf(file="./Documents/EAU_CITE-seq/ic-dotplot.pdf")
DotPlot(ics, features = feats) + 
  coord_flip() + 
  RotatedAxis() +
  theme(axis.text.y = element_text(face="italic"))
dev.off()

#CD4 T cells
cd4 <- subset(ics, idents="CD4 T cells")
DefaultAssay(cd4) <- "RNA"
cd4 <- SCTransform(cd4, method="glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"))

cd4 <- RunPCA(cd4, npcs = 50)
ElbowPlot(cd4, ndims = 50)

DefaultAssay(cd4) <- "SCT"
cd4 <- RunUMAP(cd4, dims = 1:10)
cd4 <- FindNeighbors(cd4, dims = 1:10)
cd4 <- FindClusters(cd4, resolution = seq(0.3, 1.5, 0.1))

clustree(cd4, prefix="SCT_snn_res.",
         node_colour="sc3_stability", 
         node_size_range=c(4,10), 
         edge_arrow=FALSE, 
         show_axis=TRUE)


Idents(cd4) <- "SCT_snn_res.0.8"

DefaultAssay(cd4) <- "RNA"
cd4 <- NormalizeData(cd4) %>% ScaleData()


FindAllMarkers(cd4, only.pos = T, min.pct = 0.25)

pdf("./Documents/EAU_CITE-seq/Plots/cd4-umap.pdf", height=10, width=10)
DimPlot(cd4) + NoAxes()
dev.off()

pdf("./Documents/EAU_CITE-seq/Plots/cd4-violin.pdf", height=10, width=10)
FeaturePlot(cd4, features=c("Ifng", "Il17a"))
dev.off()


data$ic_cluster <- as.character(Idents(data))
data$ic_cluster[Cells(ics)] <- paste(Idents(ics))


saveRDS(data, file="integrated_data.RDS")

DimPlot(data, group.by = "ic_cluster")

