###### maenli scRNAseq x 9.5 10.5 limb #####

library(DropletUtils)
library(Seurat)
library(Matrix)
library(dplyr)
library(biomaRt)

setwd("D:/Allou/")

## gene id info 

geneinfo <- read.table("GSE149368_genes.txt", header=T) # this list has the version id after the point

geneinfo$geneID <- gsub('\\..+$', '', geneinfo$geneID) # remove version from the list

mart <- useDataset("mmusculus_gene_ensembl", useEnsembl("ensembl", version = 98)) # version 98 if not there will be missing IDs
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "mgi_symbol"),
  values=geneinfo$geneID,
  mart=mart)

allou <- readMM("GSE149368_umi_count_matrix.mtx")

barcodes <- read.table("GSE149368_barcode_metadata.txt", header = T)

write10xCounts("D:/Allou/out/", allou,barcodes = barcodes$barcode,
               gene.symbol=genes[,2],gene.id=genes[,1])


###
maenli <- CreateSeuratObject(counts = allou,
                             project = "95vs105",assay = "RNA",
                             meta.data = "GSE149368_barcode_metadata.txt", min.cells = 3, min.features = 2000)



