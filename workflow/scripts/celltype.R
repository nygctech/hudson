library(Seurat)
library(tidyverse)
library(argparse)
parser <- ArgumentParser(description= '')
parser$add_argument('--input', '-i', help= 'Input file')
parser$add_argument('--output', '-o', help= 'Output file')
xargs <- parser$parse_args()

prot <- read_csv(xargs$input) %>% 
  column_to_rownames("...1")

so <- CreateSeuratObject(t(prot), assay = "prot")
so <- NormalizeData(so, normalization.method = "CLR")
so <- ScaleData(so) %>% RunPCA(features = rownames(so))
so <- FindNeighbors(so, dims = 1:8)
so <- RunUMAP(so, dims = 1:8)
so <- FindClusters(so, dims = 1:8, resolution = 0.2, verbose = FALSE)
write_csv(so@meta.data %>% rownames_to_column("id"), xargs$output)
