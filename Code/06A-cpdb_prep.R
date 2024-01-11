library(biomaRt)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(SeuratDisk)

#Load Seurat object
data <- readRDS("path/to/integrated_data.RDS")

# Preparing genes for CellPhoneDB as human gene ID input is required
hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get gene names in dataset
mgenes <- rownames(data)

# Get gene info from biomaRt
mgene_info <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "mgi_id"), filters = "mgi_symbol", mart = mm, values = mgenes)

# JAX table of orthologues for looking up human gene names
lookup <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")

# Filter lookup table for genes in dataset
lookup_filtered <- lookup[lookup$Mouse.MGI.ID %in% mgene_info$mgi_id,]

# Rename column for merging
colnames(lookup_filtered)[6] <- "mgi_id"

# merge to add DB.Class.Key into mgene_info
mgene_info <- merge(mgene_info, lookup_filtered[,c(1,4,6)], by="mgi_id")

# get DB.Class.Key
keys <- lookup_filtered[order(match(lookup_filtered$Mouse.MGI.ID, mgene_info$mgi_id)),1]
hlook <- filter(lookup, lookup$Common.Organism.Name == "human")


orthos <- merge(mgene_info, hlook[,c(1,4,7)], by="DB.Class.Key")
colnames(orthos)

# Function for converting gene names - if multiple orthologues, use uppercase version if present, first entry if not
z <- lapply(mgenes, function(x){
  ifelse(nrow(orthos[orthos$mgi_symbol==x,]) < 1,
         yes=toupper(x),
         no=ifelse(nrow(orthos[orthos$mgi_symbol==x,]) == 1,
                   yes=orthos[orthos$mgi_symbol==x,"Symbol.y"],
                   no=ifelse(toupper(x) %in% orthos[orthos$mgi_symbol==x,"Symbol.y"],
                             yes=toupper(x),
                             no=orthos[orthos$mgi_symbol==x,"Symbol.y"][1])))
})

z <- unlist(z)

# Create dataframe with mouse gene names alongside human orthologue
mg_list <- data.frame(Mouse.gene.name=mgenes, Human.gene.name=z)

# Rename Muller glia to avoid Unicode errors with CPDB
muller <- subset(data, ic_cluster == "Müller Glia")
muller$ic_cluster <- "Muller Glia"

data$ic_cluster[Cells(muller)] <- muller$ic_cluster

data$ic_cluster

samples <- SplitObject(data, split.by = "orig.ident")


# Extracting files for CellPhoneDB input
h_data <- samples[[1]]@assays$RNA@data
e1_data <- samples[[2]]@assays$RNA@data
e2_data <- samples[[3]]@assays$RNA@data

rownames(h_data) <- mg_list[,2]
rownames(e1_data) <- mg_list[,2]
rownames(e2_data) <- mg_list[,2]

h_meta <- cbind(rownames(samples[[1]]@meta.data), samples[[1]]@meta.data[, "cell_ids", drop=F])
e1_meta <- cbind(rownames(samples[[2]]@meta.data), samples[[2]]@meta.data[, "ic_cluster", drop=F])
e2_meta <- cbind(rownames(samples[[3]]@meta.data), samples[[3]]@meta.data[, "ic_cluster", drop=F])


healthy_conv <- CreateSeuratObject(counts=h_data, meta.data = samples[[1]]@meta.data, row.names = mg_list[,2])
eau1_conv <- CreateSeuratObject(counts=e1_data, meta.data = samples[[2]]@meta.data, row.names = mg_list[,2])
eau2_conv <- CreateSeuratObject(counts=e2_data, meta.data = samples[[3]]@meta.data, row.names = mg_list[,2])

SaveH5Seurat(healthy_conv, paste(wdir, "CPDB/", "healthy", sep=""))
SaveH5Seurat(eau1_conv, paste(wdir, "CPDB/", "eau1", sep=""))
SaveH5Seurat(eau2_conv, paste(wdir, "CPDB/", "eau2", sep=""))

Convert(paste(wdir, "CPDB/", "healthy.h5seurat", sep=""), dest = "h5ad")
Convert(paste(wdir, "CPDB/", "eau1.h5seurat", sep=""), dest = "h5ad")
Convert(paste(wdir, "CPDB/", "eau2.h5seurat", sep=""), dest = "h5ad")

write.table(h_meta, paste(wdir, "CPDB/", "healthy_metadata.txt", sep=""), sep = "\t", quote=F, row.names = F)
write.table(e1_meta, paste(wdir, "CPDB/", "eau1_metadata.txt", sep=""), sep = "\t", quote=F, row.names = F)
write.table(e2_meta, paste(wdir, "CPDB/", "eau2_metadata.txt", sep=""), sep = "\t", quote=F, row.names = F)

