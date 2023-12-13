suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(magrittr)
    library(SingleCellExperiment)
    library(scater)
    library(flexmix)
    library(splines)
    library(BiocParallel)
    library(biomaRt)
    library(miQC)
    library(Seurat)
    library(SeuratDisk)
})
mtx <- Read10X(data.dir = "scn/")
## Initialize the Seurat object with the raw (non-normalized data).
srt <- CreateSeuratObject(counts = mtx,
                                   project = "P60_SCN",
                                   min.cells = 0,
                                   min.features = 200)
srt$age <- "P60"
srt$study_id <- "wen_2020"
srt$sex <- NA
srt$tech <- "10xv2"

sce <- as.SingleCellExperiment(srt)
# mouse_mart <- useEnsembl("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
# id_map <- getBM(values = rownames(sce),
#                 filters = "ensembl_gene_id",
#                 attributes = c("ensembl_gene_id", "chromosome_name"),
#                 mart = mouse_mart)
# mt_genes <- subset(id_map,id_map$chromosome_name=="MT")$ensembl_gene_id
# feature_ctrls <- list(mito = mt_genes)

mt_genes <- grepl("^mt-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])

feature_ctrls

sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

plotMetrics(sce)

model <- mixtureModel(sce)
plotModel(sce, model)
plotFiltering(sce, model)
model2 <- mixtureModel(sce, model_type = "spline")
plotModel(sce, model2)
plotFiltering(sce, model2)
plotFiltering(sce, model2, posterior_cutoff = 0.9)
sce <- filterCells(sce, model2, posterior_cutoff = 0.9)

srt <- as.Seurat(sce)
glimpse(srt@meta.data)

SaveH5Seurat(srt, filename = "wen2020_scn.h5Seurat")
Convert("wen2020_scn.h5Seurat", dest = "h5ad")

