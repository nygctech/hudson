library(Seurat)
library(hdf5r)
library(tidyverse)
library(argparse)
library(sva)
# parse
parser <- ArgumentParser(description= 'cluster and assign cell types')
parser$add_argument('--input_h5', '-i', help = 'Input H5AD file containing protein intensity')
parser$add_argument('--output', '-o', help = 'Output file')
parser$add_argument('--ref', '-r', help = 'reference object')
parser$add_argument('--typecol', '-c', help = 'celltype column in reference metadata')
parser$add_argument('--genedict', '-d', help = 'gene dictionary if protein/antibody name is different in reference')
xargs <- parser$parse_args()
# load
file_h5 <- H5File$new(xargs$input_h5, mode = "r+")
mat <- file_h5[["X"]][,]
colnames(mat) <- file_h5[["obs/_index"]][]
rownames(mat) <- file_h5[["var/_index"]][]
file_h5$close_all()

# convert names
prot <- mat[str_detect(rownames(mat), "^[A-Z]"),]
if (!is.null(xargs$genedict)) {
  genedict <- read_csv(xargs$genedict) %>% column_to_rownames("prot")
  rownames(prot) <- sapply(rownames(prot), function(x) {
    if (x %in% rownames(genedict)) {
      return(genedict[x, "ref"])
    } else {
      return(x)
    }
  })
}

#rownames(prot) <- rownames(prot) %>% str_to_title() %>%
#  ifelse(. == "Lmn1b", "Lmnb1", .) %>%
#  ifelse(. == "Iba1", "Aif1", .) %>%
#  ifelse(. == "Nfh", "Nefh", .)

param <- mat[!str_detect(rownames(mat), "^[A-Z]"),]
# combat
#mod <- model.matrix(~as.factor(label), data = as.data.frame(t(param)))
#mod0 <- model.matrix(~1,data = as.data.frame(t(param)))
#svobj <- sva(prot, mod, mod0, n.sv = 1, method = "two-step")
#if (svobj$n.sv > 0) {
#  prot2 <- ComBat(prot, batch = svobj$sv)
#} else {
  prot2 <- prot
#}
#svobj <- sva(param, mod, mod0, n.sv = 1, method = "two-step")
#if (svobj$n.sv > 0) {
#  param2 <- ComBat(param, batch = svobj$sv)
#} else {
  param2 <- param
#}
# prot
so <- CreateSeuratObject(prot2, assay = "prot")
so <- NormalizeData(so, normalization.method = "CLR")
so@assays$prot@data <- asinh(prot / 5)
so <- ScaleData(so) %>% RunPCA(features = rownames(so))
# dr
so <- FindNeighbors(so, dims = 1:ncol(so@reductions$pca@cell.embeddings))
so <- RunUMAP(so, dims = 1:ncol(so@reductions$pca@cell.embeddings))
so <- FindClusters(so, dims = 1:ncol(so@reductions$pca@cell.embeddings), resolution = 0.1, verbose = FALSE)
# pheno
res <- CytoTree::Rphenograph(t(so@assays$prot@data))
so$phenograph <- res[[2]]$membership
# param
so[["param"]] <- CreateAssayObject(param2)
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

# write h5ad
file_h5 <- H5File$new(xargs$input_h5, mode = "r+")
dtype_string_utf8 <- H5T_STRING$new()$set_size(Inf)$set_cset("UTF-8")
#file_h5[["obs"]][["celltype"]] <- as.character(so$type)
file_h5[["obs"]]$create_dataset(so$type, name = "celltype", dtype = dtype_string_utf8)
h5attr(file_h5[["obs"]], "column-order") <- 
  tryCatch(c(h5attr(file_h5[["obs"]], "column-order"), c("celltype")) %>% unique(),
           error=function(e) c("celltype"))
file_h5$close_all()

}
