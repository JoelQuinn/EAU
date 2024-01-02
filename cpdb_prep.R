library(biomaRt)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)

#Load Seurat object
data <- readRDS("Documents/EAU_CITE-seq/integrated_data.RDS")

# Preparing genes for CellPhoneDB as requires human gene ID input ----
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

mg_list <- data.frame(Mouse.gene.name=mgenes, Human.gene.name=z)


samples <- SplitObject(data, split.by = "orig.ident")

# Extracting files for CellPhoneDB input
h_data <- samples[[1]]@assays$RNA@data
e1_data <- samples[[2]]@assays$RNA@data
e2_data <- samples[[3]]@assays$RNA@data

rownames(h_data) <- mg_list[,2]
rownames(e1_data) <- mg_list[,2]
rownames(e2_data) <- mg_list[,2]

write.table(h_data, "healthy_counts.txt", sep="\t", quote=F)
write.table(e1_data, "eau1_counts.txt", sep="\t", quote=F, col.names = colnames(e1_data))
write.table(e2_data, "eau2_counts.txt", sep="\t", quote=F, col.names = colnames(e2_data))

colnames(e1_data)

samples[[2]]@meta.data

h_data

h_meta <- cbind(rownames(samples[[1]]@meta.data), samples[[1]]@meta.data[, "cell_ids", drop=F])
e1_meta <- cbind(rownames(samples[[2]]@meta.data), samples[[2]]@meta.data[, "ic_cluster", drop=F])
e2_meta <- cbind(rownames(samples[[3]]@meta.data), samples[[3]]@meta.data[, "ic_cluster", drop=F])

write.table(h_meta, "healthy_metadata.txt", sep = "\t", quote=F, row.names = F)
write.table(e1_meta, "eau1_metadata.txt", sep = "\t", quote=F, row.names = F)
write.table(e2_meta, "eau2_metadata.txt", sep = "\t", quote=F, row.names = F)