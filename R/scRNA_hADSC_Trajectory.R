
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

`%ni%` <- function(a, b) return(! a %in% b)

# WOT ---------------------------------------------------------------------

allSo <- readRDS("middata/hADSC_0618_scRNA/allSo_qc.rds")

# prepare input

usedSample <- c(
  "ADSC",
  "S1D4", "S1D8",
  "S2D4", "S2D8", "S2D12", "S2D16"
)

sampleTime <- c(
  0, 4, 8,
  8 + c(4, 8, 12, 16),
)

names(sampleTime) <- usedSample

allSo <- allSo[, allSo$sample %in% usedSample]

allSo$time <- sampleTime[allSo$sample]
table(allSo@meta.data[, c("time", "sample")])

wotExp <- allSo[["RNA"]]@data[VariableFeatures(allSo), ] %>% as.matrix %>% t
wotExp <- cbind(rownames(wotExp), wotExp)
colnames(wotExp)[1] <- "id"

write.table(wotExp, "middata/hADSC_scRNA/wot/wotExp.txt", sep = "\t", row.names = F, quote = F)

wotDay <- as.matrix(allSo$time)
rownames(wotDay) <- colnames(allSo)
wotDay <- cbind(rownames(wotDay), wotDay)
colnames(wotDay) <- c("id", "day")
write.table(wotDay, "middata/hADSC_scRNA/wot/wotDay.txt", sep = "\t", quote = F, row.names = F)

# perform WOT analysis according to developer instructions (https://broadinstitute.github.io/wot/)

# identify S2D16 celltype

end <- allSo[, allSo$sample == "S2D16"]
end %<>% FindVariableFeatures
end %<>% ScaleData %>% RunPCA
end %>% ElbowPlot(ndims = 50)
end %<>% FindNeighbors(dims = 1:20)
end %<>% FindClusters(resolution = 1)
end %<>% RunUMAP(dims = 1:20)

DimPlot(end, group.by = "seurat_clusters", label = T)

FeaturePlot(end, c("SALL4", "LIN28A", "MSX1", "MSX2", "HOXB9"))
FeaturePlot(end, c("COL1A2", "DCN"))

end$celltype <- "other"
end$celltype[end$seurat_clusters %in% c()] <- "IPC"

# the celltypes were used to calculate right probability

allSo$IPC <- readRDS("middata/hADSC_scRNA/IPC_score.rds")[colnames(allSo)]

saveRDS(allSo, "middata/hADSC_scRNA/allSo_traj.rds")

# identify right trajectory -----------------------------------------------

allSo <- readRDS("middata/hADSC_scRNA/allSo_traj.rds")

allSo$traj <- "F"
allSo$traj[allSo$sample == "prime"] <- "prime"
allSo$traj[allSo$sample == "H1"] <- "ES"

# define right trajectory according to IPC score
allSo$traj[allSo$seurat_clusters %in% c()] <- "R0"
allSo$traj[allSo$seurat_clusters %in% c()] <- "R1"
allSo$traj[allSo$seurat_clusters %in% c()] <- "R2"
allSo$traj[allSo$seurat_clusters %in% c()] <- "R3"
allSo$traj[allSo$seurat_clusters %in% c()] <- "R4"

saveRDS(allSo, "middata/hADSC_scRNA/allSo_traj.rds")

# integrate scRNA and scATAC+RNA ------------------------------------------

# after scATAC+RNA_hADSC_PreProcess.R

sc <- readRDS("middata/hADSC_scRNA/allSo_traj.rds")
mo <- readRDS("middata/hADSC_scMultiome/allSo_qc.rds")

usedGene <- intersect(rownames(sc), rownames(mo))

countData <- cbind(
  sc[["RNA"]]@counts[usedGene, ],
  mo[["RNA"]]@counts[usedGene, ]
)
colnames(countData) <- c(
  str_c("sc_", colnames(sc)),
  str_c("mo_", colnames(mo))
)

scmo <- CreateSeuratObject(countData)

scmo$sample <- c(
  sc$sample %>% as.vector(),
  mo$sample %>% as.vector()
)
scmo$celltype <- c(
  sc$traj %>% as.vector(),
  mo$sample %>% as.vector()
)
scmo$data <- c(
  rep("sc", ncol(sc)),
  rep("mo", ncol(mo))
)

scmo %<>% NormalizeData()

scmo %<>% FindVariableFeatures()
scmo %<>% ScaleData() %>% RunPCA()
scmo %>% ElbowPlot(ndims = 50)

sc <- scmo[, scmo$data == "tri"]
mo <- scmo[, scmo$data == "mo"]

library(symphony)

ref <- buildReference(
  sc[["RNA"]]@data,
  sc@meta.data,
  vars = c("Phase"),
  K = 100,
  verbose = T,
  do_umap = T,
  do_normalize = F,
  vargenes_groups = NULL,
  topn = 1000,
  save_uwot_path = "sc_umap.model"
)
saveRDS(ref, "middata/hADSC_scRNA/sc_symphonyRef.rds")

sc[["symphony"]] <- CreateDimReducObject(embeddings = ref$Z_corr %>% t, assay = "RNA")
sc[["umap"]] <- CreateDimReducObject(
  embeddings = ref$umap$embedding %>% set_rownames(colnames(sc)) %>% set_colnames(c("UMAP_1", "UMAP_2")),
  assay = "RNA")

moProj <- mapQuery(
  exp_query = mo[["RNA"]]@data,
  metadata_query = mo@meta.data,
  ref_obj = ref,
  vars = c("Phase"),
  do_normalize = F,
  do_umap = T,
  sigma = 0.1
)

mo[["symphony"]] <- CreateDimReducObject(embeddings = moProj$Z %>% t, assay = "RNA")
mo[["umap"]] <- CreateDimReducObject(
  embeddings = moProj$umap %>% set_rownames(colnames(mo)) %>% set_colnames(c("UMAP_1", "UMAP_2")),
  assay = "RNA")

scmo[["symphony"]] <- CreateDimReducObject(rbind(sc[["symphony"]]@cell.embeddings, mo[["symphony"]]@cell.embeddings)[colnames(scmo), ], assay = "RNA")
scmo[["umap"]] <- CreateDimReducObject(rbind(sc[["umap"]]@cell.embeddings, mo[["umap"]]@cell.embeddings)[colnames(trimo), ], assay = "RNA")

plotMeta(scmo, "sample")
plotMeta(scmo, "data")

scmo %<>% FindNeighbors(dims = 1:20, reduction = "symphony")
scmo %<>% FindClusters(resolution = 1)

saveRDS(scmo, "middata/hADSC_scRNA/scmo.rds")

# annotate sc to mo -------------------------------------------------------

scPCA <- scmo[["symphony"]]@cell.embeddings[scmo$data == "sc", 1:20]
moPCA <- scmo[["symphony"]]@cell.embeddings[scmo$data == "mo", 1:20]

mat1 <- scPCA %*% t(moPCA)

scL <- apply(scPCA, 1, crossprod)
moL <- apply(moPCA, 1, crossprod)

mat2 <- matrix(triL, nrow = nrow(mat1), ncol = ncol(mat1))
mat3 <- matrix(moL, nrow = ncol(mat1), ncol = nrow(mat1)) %>% t

distMat <- mat2 + mat3 - 2*mat1
remove(mat1, mat2, mat3); gc()
saveRDS(distMat, "middata/hADSC_scRNA/mo_sc_distMat.rds")
distMat <- readRDS("middata/hADSC_scRNA/mo_sc_distMat.rds")

nnName <- apply(distMat, 2, function(x) {names(sort(x)[1:20])})
saveRDS(nnName, "middata/hADSC_scRNA/mo_nnName.rds")
nnName <- readRDS("middata/hADSC_scRNA/mo_nnName.rds")

nnDist <- apply(distMat, 2, function(x) {sort(x)[1:20]})
saveRDS(nnDist, "middata/hADSC_scRNA/mo_nnDist.rds")

mat1 <- tcrossprod(scPCA, scPCA)
mat2 <- matrix(triL, nrow = nrow(mat1), ncol = ncol(mat1))

distMat <- mat2 + t(mat2) - 2*mat1
remove(mat1, mat2); gc()
saveRDS(distMat, "middata/hADSC_scRNA/sc_self_distMat.rds")

sc %<>% FindNeighbors(dims = 1:20, reduction = "symphony")
snn <- summary(sc@graphs$RNA_snn) %>% as.data.table()

min(table(snn$i))


library(progress)
pb <- progress_bar$new(
  format = "running [:bar]:percent | eta :eta | already :elapsed | i = :current",
  total = ncol(tri),
  clear = F,
  stream = stdout())

binWid <- c()
for(p in 1:ncol(sc)) {
  pb$tick()
  neiber <- snn[i == p][order(x)]
  
  if(nrow(neiber) <= 20) {
    binWid <- c(binWid, distMat[p, neiber[, j]] %>% sqrt %>% mean)
  } else if(neiber[20, x] == neiber[21, x]) {
    n <- sum(neiber$x <= neiber[20, x])
    ndrop <- n - 20
    id <- neiber[x == neiber[20, x], j]
    idDist <- distMat[p, id] %>% sort()
    binWid <- c(binWid, c(idDist[-(1:ndrop)], distMat[p, neiber[x < neiber[20, x], j]], idDist) %>% sqrt %>% mean)
  } else {
    binWid <- c(binWid, distMat[p, neiber[1:20, j]] %>% sqrt %>% mean)
  }
}
saveRDS(binWid, "middata/hADSC_scRNA/sc_binWid.rds")
binWid <- readRDS("middata/hADSC_scRNA/sc_binWid.rds")

nnName[, 1:5]
moBinWid <- apply(nnName, 2, function(x) mean(binWid[match(x, colnames(tri))]))
saveRDS(moBinWid, "middata/hADSC_scRNA/mo_binWid.rds")

nnDist[, 1:5]
nnW <- sweep(-nnDist/2, 2, (moBinWid)^2, "/")
nnW[, 1:10]
nnW <- apply(nnW, 2, function(x) {w <- exp(x); w/sum(w)})
saveRDS(nnW, "middata/hADSC_scRNA/mo_nnW.rds")
nnW <- readRDS("middata/hADSC_scRNA/mo_nnW.rds")

for(i in str_c("R", 1:4)) {
  message(i)
  
  tri@meta.data[, str_c("is", i)] <- ifelse(tri$celltype == i, 1, 0)
  nnFeat <- tri@meta.data[, str_c("is", i)][match(nnName %>% as.vector, colnames(tri))]
  
  mo@meta.data[, str_c("is", i)] <- matrix(nnFeat * as.vector(nnW), nrow = 20) %>% colSums()
  mo@meta.data[, str_c("is", i)] %<>% {.[is.na(.)] <- 0; .}
}

mo$traj <- "F"
mo$traj[mo$sample == "ADSC"] <- "R0"
mo$traj[mo$sample != "ADSC" & mo$isR1 > 0.5] <- "R1"
mo$traj[mo$sample != "ADSC" & mo$isR2 > 0.5] <- "R2"
mo$traj[mo$sample != "ADSC" & mo$isR3 > 0.5] <- "R3"
mo$traj[mo$sample != "ADSC" & mo$isR4 > 0.5] <- "R4"

scmo$celltype <- c(
  sc$traj %>% as.vector(),
  mo$traj %>% as.vector()
)

saveRDS(mo, "middata/hADSC_scMultiome/allSo_traj.rds")
saveRDS(scmo, "middata/hADSC_scRNA/scmo.rds")

# PAGA --------------------------------------------------------------------

cellMeta <- scmo@meta.data[, c("part", "sample", "traj")]
geneMeta <- GetAssay(scmo)[[]]
usedGene <- VariableFeatures(scmo)

library(reticulate)
sc <- import("scanpy")
adata <- sc$AnnData(
  X = scmo[["RNA"]]@data[usedGene, ] %>% Matrix::t(),
  obs = cellMeta,
  var = geneMeta[usedGene, ]
)
adata$obsm$update(X_pca = Embeddings(scmo, "lsi"))
adata$obsm$update(X_umap = Embeddings(scmo, "umap"))

adata$write_h5ad("middata/hADSC_scRNA/allAnn.h5ad")

# perform PAGA analysis (https://github.com/theislab/paga) to get FA coordinates of trajectory

adata <- sc$read_h5ad("middata/hADSC_scRNA/allAnn_paga.h5ad")

fa_emb <- adata$obsm['X_draw_graph_fa']
colnames(fa_emb) <- c("FA_1", "FA_2")
row.names(fa_emb) <- adata$obs_names$to_list()

scmo[["fa"]] <- CreateDimReducObject(embeddings = fa_emb, key = "FA_", assay = "RNA")

DimPlot(scmo, group.by = "sample", reduction = "fa", label = T)

saveRDS(scmo, "middata/hADSC_scRNA/scmo.rds")

