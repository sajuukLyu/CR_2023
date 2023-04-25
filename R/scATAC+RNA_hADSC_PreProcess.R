
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu, Cheng Li Lab, Peking University
# email: lvyulin@pku.edu.cn
# date: 2023.4.24
#
# ---

# load package ------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(data.table)

library(Seurat)
library(SingleCellExperiment)
library(scran)
library(ArchR)

source("R/utils.R")

`%ni%` <- function(a, b) return(! a %in% b)

# load data ---------------------------------------------------------------

usedSample <- c(
  "ADSC", "S1D9", "S2D8", "S2D16", "S3D10"
)

dataList <- map(str_c("data/hADSC_scMultiome/", usedSample, "/filtered_feature_bc_matrix.h5"), Read10X_h5)
names(dataList) <- usedSample

dataList <- map(dataList, ~ .x$`Gene Expression`)
dataList <- imap(dataList, ~ {set_colnames(.x, str_c(.y, "_", colnames(.x)))})
map_int(dataList, ncol)

countData <- do.call(cbind, dataList)
allSo <- CreateSeuratObject(countData)

saveRDS(allSo, "middata/hADSC_scMultiome/allSo_raw.rds")

# scRNA quality control ---------------------------------------------------

allSo <- readRDS("middata/hADSC_scMultiome/allSo_raw.rds")
allSo$sample <- str_extract(colnames(allSo), ".*_") %>% str_remove("_$")
table(allSo$sample)

mtGene <- rownames(allSo) %>% str_subset("^MT-")
allSo$mt <- PercentageFeatureSet(allSo, features = mtGene)

allSo <- allSo[, allSo$nFeature_RNA > 500 & allSo$nCount_RNA > 1000 & allSo$mt < 20]
allSo %<>% NormalizeData

# check doublet
soList <- SplitObject(allSo, "sample")
soList <- map(soList, FindVariableFeatures)
SCElist <- map(soList, ~ SingleCellExperiment(assays = list(counts = .x[["RNA"]]@counts, logcounts = .x[["RNA"]]@data)))
varGeneList <- map(soList, VariableFeatures)
dblDen <- map2(SCElist, varGeneList, ~ doubletCells(.x, subset.row = intersect(rownames(.x), .y), d = 20))
dblDen <- map(dblDen, ~ log10(.x + 1))
allSo$dbl <- unlist(dblDen) %>% unname

allSo %<>% FindVariableFeatures
allSo %<>% ScaleData %>% RunPCA
allSo %>% ElbowPlot(ndims = 50)
allSo %<>% FindNeighbors(dims = 1:20)
allSo %<>% FindClusters(resolution = 1)
allSo %<>% RunUMAP(dims = 1:20)

DimPlot(allSo, group.by = "sample", label = T)

VlnPlot(allSo, "nCount_RNA", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(allSo, "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(allSo, "dbl", group.by = "seurat_clusters", pt.size = 0)

badClu <- c() # low quality clusters
allSo <- allSo[, allSo$seurat_clusters %ni% badClu]

saveRDS(allSo, "middata/hADSC_scMultiome/allSo_qc.rds")

# scATAC quality control --------------------------------------------------

addArchRGenome("hg19")
addArchRThreads(threads = 12)

usedSample <- c("ADSC", "S1D9", "S2D8", "S2D16", "S3D10")
inputFiles <- str_c("data/hADSC_scMultiome/", usedSample, "/atac_fragments.tsv.gz")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = usedSample,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = T,
  addGeneScoreMat = T,
  force = T
)

ArrowFiles <- list.files(".", "arrow")

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

ao <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "middata/hADSC_scMultiome/atac",
  copyArrows = T
)

# perform QC
ao <- addIterativeLSI(ao, force = T)
ao <- addUMAP(
  ArchRProj = ao,
  nNeighbors = 30,
  minDist = 0.5,
  force = T
)
ao <- addClusters(
  input = ao,
  resolution = 1,
  force = T
)

plotATACmeta(ao, "Sample")
plotATACmeta(ao, "Clusters")
table(ao$Clusters)

plotATACmeta(ao, "TSSEnrichment", "c", lb_fs = "Clusters")
plotATACmeta(ao, "nFrags", "c", lb_fs = "Clusters")
plotATACmeta(ao, "DoubletScore", "c", lb_fs = "Clusters")
plotATACmeta(ao, "DoubletEnrichment", "c", lb_fs = "Clusters")

getCellColData(ao, c("Clusters", "nFrags")) %>% as.data.table %>% ggplot(aes(x = Clusters, y = nFrags)) +
  geom_violin(scale = "width") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
getCellColData(ao, c("Clusters", "TSSEnrichment")) %>% as.data.table %>% ggplot(aes(x = Clusters, y = TSSEnrichment)) +
  geom_violin(scale = "width") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
getCellColData(ao, c("Clusters", "DoubletEnrichment")) %>% as.data.table %>% ggplot(aes(x = Clusters, y = DoubletEnrichment)) +
  geom_violin(scale = "width") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

badClu <- str_c("C", c()) # low quality clusters
ao <- ao[ao$Clusters %ni% badClu, ]

ao <- saveArchRProject(ao, outputDirectory = "middata/hADSC_scMultiome/atac", load = T)

# filter scRNA+scATAC -----------------------------------------------------

allSo <- readRDS("middata/hADSC_scMultiome/allSo_qc.rds")
ao <- loadArchRProject("middata/hADSC_scMultiome/atac")

# keep cells passed both scRNA and scATAC quality control
atacName <- str_replace(ao$cellNames, "#", "_")
allSo$hasATAC <- colnames(allSo) %in% atacName

allSo <- allSo[, allSo$hasATAC]
ao <- ao[match(colnames(allSo), atacName), ]

saveRDS(allSo, "middata/hADSC_scMultiome/allSo_qc.rds")
ao <- saveArchRProject(ao, outputDirectory = "middata/hADSC_scMultiome/atac", load = T)

