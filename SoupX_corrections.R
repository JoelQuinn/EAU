library(SoupX)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(DropletUtils)
library(clustree)

# Dir names
h_dir <- "D:/CITE-seq/REX/WTCHG_917226_TY73/outs"
res1_dir <- "D:/CITE-seq/REX/WTCHG_917226_TY74/outs"
res2_dir <- "D:/CITE-seq/REX/WTCHG_917226_TY75/outs"
eau1_dir <- "D:/CITE-seq/REX/WTCHG_917226_TY76/outs"
eau2_dir <- "D:/CITE-seq/REX/WTCHG_917226_TY77/outs"

# Untreated (healthy)
unt_raw <- Read10X(paste(h_dir, "/raw_feature_bc_matrix", sep=""))
unt_filt <- Read10X(paste(h_dir, "/filtered_feature_bc_matrix", sep=""))

# Resistant 1
res1_raw <- Read10X(paste(res1_dir, "/raw_feature_bc_matrix", sep=""))
res1_filt <- Read10X(paste(res1_dir, "/filtered_feature_bc_matrix", sep=""))

# Resistant 2
res2_raw <- Read10X(paste(res2_dir, "/raw_feature_bc_matrix", sep=""))
res2_filt <- Read10X(paste(res2_dir, "/filtered_feature_bc_matrix", sep=""))

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

# Clustree plotting function
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
wt <- CreateSeuratObject(counts = unt_filt)
wtSoup <- SoupChannel(tod = unt_raw, toc=unt_filt)

res1 <- CreateSeuratObject(counts=res1_filt)
res1soup <- SoupChannel(res1_raw, res1_filt)

res2 <- CreateSeuratObject(counts=res2_filt)
res2soup <- SoupChannel(res2_raw, res2_filt)

eau1 <- CreateSeuratObject(counts=eau1_filt)
eau1soup <- SoupChannel(eau1_raw, eau1_filt)

eau2 <- CreateSeuratObject(counts=eau2_filt)
eau2soup <- SoupChannel(eau2_raw, eau2_filt)

#----
#healthy
wt <- transPCA(wt)
ElbowPlot(wt, ndims = 75)
wt <- umap_cluster(wt, dims = 40)
quickclustree(wt)
wt_out <- soupcorrect(wt, wtSoup, 0.8)
write10xCounts(path = "SoupX_TY74", x=wt_out)


# resistant 1
res1 <- transPCA(res1)
ElbowPlot(res1, ndims = 75)
res1 <- umap_cluster(res1, dims = 40)
quickclustree(res1)
res1_out <- soupcorrect(res1, res1soup, 0.8)
write10xCounts(path = "SoupX_TY75", x=res1_out)



# Res2
res2 <- transPCA(res2)
ElbowPlot(res2, ndims = 75)
res2 <- umap_cluster(res2, dims = 35)
quickclustree(res2)
res2_out <- soupcorrect(res2, res2soup)
write10xCounts(path = "SoupX_TY76", x=res2_out)

#EAU1
eau1 <- transPCA(eau1)
ElbowPlot(eau1, ndims = 75)
eau1 <- umap_cluster(eau1, dims = 50)
quickclustree(eau1)
eau1_out <- soupcorrect(eau1, eau1soup, 1)
write10xCounts(path = "SoupX_TY77", x=eau1_out)

#EAU2
eau2 <- transPCA(eau2)
ElbowPlot(eau2, ndims = 75)
eau2 <- umap_cluster(eau2, dims = 50)
quickclustree(eau2)
eau2_out <- soupcorrect(eau2, eau2soup)
write10xCounts(path = "SoupX_TY78", x=eau2_out)
        