library(fgsea)
library(dplyr)
library(tidyr)
library(BiocGenerics)
library(msigdbr)

directory <- "path/to/DESeq2_results"

# Prepare Hallmark pathways - .rdata object downloaded from https://bioinf.wehi.edu.au/software/MSigDB/
load("Documents/EAU_CITE-seq/mouse_H_v5p2.rdata")
pathwaysH <- Mm.H

# Prepare GO:BP pathways
c5 <- msigdbr::msigdbr(species="Mus musculus", category="C5") %>%
  filter(gs_subcat == "GO:BP")

pathwaysC5 <- split(c5$entrez_gene, f=c5$gs_name)

# List of folders
folders <- list.files(directory)

# Run analysis in each folder (cell type)
# Note this can be very memory intensive for GO:BP gene sets, recommend running on HPC cluster
for (i in folders) {  
    filename <- paste(directory, "/", i, "/", i, "_ranked-genes_ENTREZ.csv", sep="")
    
    # Read file
    ranks <- read.csv(file = filename)
    # Remove NAs
    ranks <- na.omit(ranks)

    # Get gene names and rank
    gene_ranks <- ranks$gene_list
    names(gene_ranks) <- ranks$X
    
    gene_ranks <- sort(gene_ranks, decreasing = T)

    # Run fgsea on Hallmark
    fgseaResH <- fgsea(pathways = pathwaysH, stats = gene_ranks, 
                      minSize=15, maxSize=500)

    # Run fgsea on GO:BP
    fgseaResC5 <- fgsea(pathways = pathwaysC5, stats = gene_ranks, 
                        minSize=15, maxSize=500)
    
    fgseaResHOrdered <- fgseaResH[order(padj, -abs(NES)), ]
    fgseaResHOrdered <- as.data.frame(fgseaResHOrdered)
    fgseaResHOrdered <- apply(fgseaResHOrdered,2,as.character)
    
    fgseaResC5Ordered <- fgseaResC5[order(padj, -abs(NES)), ]
    fgseaResC5Ordered <- as.data.frame(fgseaResC5Ordered)
    fgseaResC5Ordered <- apply(fgseaResC5Ordered,2,as.character)

    # Create output directory
    fgsea_dir <- paste(directory, "/", i, "/", i, "_fgsea", sep="")
    dir.create(fgsea_dir)

    # Output results
    gseaH_out <- paste(fgsea_dir, "/", i, "_fgsea-HALLMARK-results-table.csv", sep="")
    gseaC5_out <- paste(fgsea_dir, "/", i, "_fgsea-C5BP-results-table.csv", sep="")
    
    write.csv(x=fgseaResHOrdered, file=gseaH_out)
    write.csv(x=fgseaResC5Ordered, file=gseaC5_out)
    
}

