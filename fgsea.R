library(fgsea)
library(dplyr)
library(tidyr)
library(BiocGenerics)
library(msigdbr)

directory <- commandArgs(trailingOnly = TRUE)

load("Documents/EAU_CITE-seq/mouse_H_v5p2.rdata")
pathwaysH <- Mm.H

c5 <- msigdbr::msigdbr(species="Mus musculus", category="C5") %>%
  filter(gs_subcat == "GO:BP")

pathwaysC5 <- split(c5$entrez_gene, f=c5$gs_name)



folders <- list.files(directory)
  
for (i in folders) {
    filename <- paste(directory, "/", i, "/", i, "_ranked-genes_ENTREZ.csv", sep="")
    ranks <- read.csv(file = filename)
    ranks <- na.omit(ranks)
    
    gene_ranks <- ranks$gene_list
    names(gene_ranks) <- ranks$X
    
    gene_ranks <- sort(gene_ranks, decreasing = T)
    
    fgseaResH <- fgsea(pathways = pathwaysH, stats = gene_ranks, 
                      minSize=15, maxSize=500)
    
    fgseaResC5 <- fgsea(pathways = pathwaysC5, stats = gene_ranks, 
                        minSize=15, maxSize=500)
    
    fgseaResHOrdered <- fgseaResH[order(padj, -abs(NES)), ]
    fgseaResHOrdered <- as.data.frame(fgseaResHOrdered)
    fgseaResHOrdered <- apply(fgseaResHOrdered,2,as.character)
    
    fgseaResC5Ordered <- fgseaResC5[order(padj, -abs(NES)), ]
    fgseaResC5Ordered <- as.data.frame(fgseaResC5Ordered)
    fgseaResC5Ordered <- apply(fgseaResC5Ordered,2,as.character)
    
    fgsea_dir <- paste(directory, "/", i, "/", i, "_fgsea", sep="")
    dir.create(fgsea_dir)
    
    gseaH_out <- paste(fgsea_dir, "/", i, "_fgsea-HALLMARK-results-table.csv", sep="")
    gseaC5_out <- paste(fgsea_dir, "/", i, "_fgsea-C5BP-results-table.csv", sep="")
    
    write.csv(x=fgseaResHOrdered, file=gseaH_out)
    write.csv(x=fgseaResC5Ordered, file=gseaC5_out)
    
}

