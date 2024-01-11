library(SoupX)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(DropletUtils)
library(clustree)

# Dir names
h_dir <- "path/to/10X_files/healthy"
eau1_dir <- "path/to/10X_files/eau1"
eau2_dir <- "path/to/10X_files/eau2"

# Healthy
h_raw <- Read10X(paste(h_dir, "/raw_feature_bc_matrix", sep=""))
h_filt <- Read10X(paste(h_dir, "/filtered_feature_bc_matrix", sep=""))

# EAU 1
eau1_raw <- Read10X(paste(eau1_dir, "/raw_feature_bc_matrix", sep=""))
eau1_filt <- Read10X(paste(eau1_dir, "/filtered_feature_bc_matrix", sep=""))

# EAU 2
eau2_raw <- Read10X(paste(eau2_dir, "/raw_feature_bc_matrix", sep=""))
eau2_filt <- Read10X(paste(eau2_dir, "/filtered_feature_bc_matrix", sep=""))

# Functions for quick processing of data for SoupX input
# SCTransform and PCA
transPCA <- function(srat, n=75){
  srat <- SCTransform(srat)
  srat <- RunPCA(srat, npcs = n)
  return(srat)
}
# UMAP dim reduction and clustering
umap_cluster <- function(srat, dims){
  srat <- RunUMAP(srat, dims=1:dims) %>%
    FindNeighbors(dims=1:dims) %>% 
    FindClusters(resolution = seq(0.4, 1.5, 0.1))
  return(srat)
}

# Clustree plotting function to choose optimal clustering resolution
quickclustree <- function(srat) {
  clustree(srat, prefix="SCT_snn_res.",
           node_colour="sc3_stability", 
           node_size_range=c(4,10), 
           edge_arrow=FALSE, 
           show_axis=TRUE)
}

# Function for running SoupX steps
soupcorrect <- function(srat, soup, cluster_resolution=0.8){
  
  mcolumn <- paste("SCT_snn_res.", cluster_resolution, sep="")
  
  srat@meta.data$seurat_clusters <- srat@meta.data[[mcolumn]]
  
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  
  soup <- setClusters(soup, setNames(meta$seurat_clusters, rownames(soup)))
  soup <- setDR(soup, umap)
  
  soup <- autoEstCont(soup)
  head(soup$soupProfile[order(soup$soupProfile$est, decreasing = T), ], n = 20)
  
  out <- adjustCounts(soup, roundToInt = TRUE)
  
  return(out)
}



#----
# loading objects
h <- CreateSeuratObject(counts = h_filt)
hsoup <- SoupChannel(tod = unt_raw, toc=h_filt)

eau1 <- CreateSeuratObject(counts=eau1_filt)
eau1soup <- SoupChannel(eau1_raw, eau1_filt)

eau2 <- CreateSeuratObject(counts=eau2_filt)
eau2soup <- SoupChannel(eau2_raw, eau2_filt)

#----
#healthy
h <- transPCA(h)
ElbowPlot(h, ndims = 75)
h <- umap_cluster(h, dims = 40)
quickclustree(h)
h_out <- soupcorrect(h, hsoup, 0.8)
write10xCounts(path = "SoupX_Healthy", x=h_out)

#EAU1
eau1 <- transPCA(eau1)
ElbowPlot(eau1, ndims = 75)
eau1 <- umap_cluster(eau1, dims = 50)
quickclustree(eau1)
eau1_out <- soupcorrect(eau1, eau1soup, 1)
write10xCounts(path = "SoupX_EAU1", x=eau1_out)

#EAU2
eau2 <- transPCA(eau2)
ElbowPlot(eau2, ndims = 75)
eau2 <- umap_cluster(eau2, dims = 50)
quickclustree(eau2)
eau2_out <- soupcorrect(eau2, eau2soup)
write10xCounts(path = "SoupX_EAU2", x=eau2_out)
        
