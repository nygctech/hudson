library(Seurat)
library(hdf5r)
library(tidyverse)
library(argparse)
# parse
parser <- ArgumentParser(description= 'cluster and assign cell types')
parser$add_argument('--input_h5', '-i', help = 'Input H5AD file containing protein intensity')
parser$add_argument('--output', '-o', help = 'Output file')
parser$add_argument('--ref', '-r', help = 'reference object')
parser$add_argument('--typecol', '-c', help = 'celltype column in reference metadata')
xargs <- parser$parse_args()
# load
file_h5 <- H5File$new(xargs$input_h5, mode = "r+")
mat <- file_h5[["X"]][,]
colnames(mat) <- file_h5[["obs/_index"]][]
rownames(mat) <- file_h5[["var/_index"]][]
prot <- mat[str_detect(rownames(mat), "_Intensity"),]
rownames(prot) <- str_remove(rownames(prot), "_Intensity")
rownames(prot) <- rownames(prot) %>% str_to_title() %>%
  ifelse(. == "Iba1", "Aif1", .) %>%
  ifelse(. == "Nfh", "Nefh", .)
param <- mat[!str_detect(rownames(mat), "_Intensity"),]
# prot
so <- CreateSeuratObject(prot, assay = "prot")
so <- NormalizeData(so, normalization.method = "CLR")
so <- ScaleData(so) %>% RunPCA(features = rownames(so))
so <- FindNeighbors(so, dims = 1:ncol(so@reductions$pca@cell.embeddings))
so <- RunUMAP(so, dims = 1:ncol(so@reductions$pca@cell.embeddings))
so <- FindClusters(so, dims = 1:ncol(so@reductions$pca@cell.embeddings), resolution = 0.1, verbose = FALSE)
# param
so[["param"]] <- CreateAssayObject(param)
DefaultAssay(so) <- "param"
so <- ScaleData(so, features = rownames(so))
so <- so %>% RunPCA(features = rownames(so), reduction.name = 'ppca')
DefaultAssay(so) <- "prot"
# celltype
if (!file.exists(xargs$ref)) {
  write_csv(so@meta.data %>% rownames_to_column("id") %>% select(id, seurat_clusters), xargs$output)
  #blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  #msg <- simpleError(blankMsg)
  #stop(msg)
} else {

ref <- readRDS(xargs$ref)
cell_type <- xargs$typecol
anchors <- suppressWarnings(FindTransferAnchors(reference = ref,
                                                query = so,
                                                features = rownames(so),
                                                reduction = "cca"))
pred <- TransferData(anchors, ref@meta.data[[cell_type]], k.weight = 10, weight.reduction = so[["pca"]], dims = 1:ncol(so@reductions$pca@cell.embeddings))
so$type <- pred$predicted.id
write_csv(so@meta.data %>% rownames_to_column("id") %>% select(id, type), xargs$output)

}
