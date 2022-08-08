library(Seurat)
library(tidyverse)
library(argparse)
# parse
parser <- ArgumentParser(description= 'cluster and assign cell types')
parser$add_argument('--input_prot', '-i', help = 'Input protein intensity')
parser$add_argument('--input_param', '-p', help = 'Input parameters')
parser$add_argument('--output', '-o', help = 'Output file')
parser$add_argument('--ref', '-r', help = 'reference object')
xargs <- parser$parse_args()
# load
prot <- read_csv(xargs$input_prot) %>% 
  column_to_rownames("...1")
colnames(prot) <- colnames(prot) %>% str_to_title() %>%
  ifelse(. == "Iba1", "Aif1", .) %>%
  ifelse(. == "Nfh", "Nefh", .)
param <- read_csv(xargs$input_param) %>%
  column_to_rownames("...1")
# prot
so <- CreateSeuratObject(t(prot), assay = "prot")
so <- NormalizeData(so, normalization.method = "CLR")
so <- ScaleData(so) %>% RunPCA(features = rownames(so))
so <- FindNeighbors(so, dims = 1:8)
so <- RunUMAP(so, dims = 1:8)
so <- FindClusters(so, dims = 1:8, resolution = 0.2, verbose = FALSE)
# param
so[["param"]] <- CreateAssayObject(t(param))
DefaultAssay(so) <- "param"
so <- ScaleData(so, features = rownames(so))
so <- so %>% RunPCA(features = rownames(so), reduction.name = 'ppca')
DefaultAssay(so) <- "prot"
# celltype
ref <- readRDS(xargs$ref)
anchors <- suppressWarnings(FindTransferAnchors(reference = ref,
                                                query = so,
                                                features = rownames(so),
                                                reduction = "cca"))
pred <- TransferData(anchors, ref$cell_type, k.weight = 10, weight.reduction = so[["pca"]], dims = 1:8)
so$type <- pred$predicted.id
write_csv(so@meta.data %>% rownames_to_column("id") %>% select(id, type), xargs$output)
