
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

# peak 2 gene -------------------------------------------------------------

allSo <- readRDS("middata/hADSC_scMultiome/allSo_qc.rds")
ao <- loadArchRProject("middata/hADSC_scMultiome/atac")

# add gene expression data to ArchR object
seRNA <- import10xFeatureMatrix(
  input = c("data/hADSC_scMultiome/ADSC/filtered_feature_bc_matrix.h5"),
  names = c("ADSC")
)

rr <- seRNA@rowRanges
seqlevels(rr) <- seqlevels(rr)[1:24]
usedGene <- intersect(rr$name, rownames(allSo))
rr <- rr[usedGene, ]

countData <- allSo[["RNA"]]@counts[usedGene, ]
colnames(countData) <- ao$cellNames

se <- SummarizedExperiment(
  assays = SimpleList(counts = countData),
  rowRanges = rr
)
ao <- addGeneExpressionMatrix(input = ao, seRNA = se, force = T)
ao <- saveArchRProject(ao, outputDirectory = "middata/hADSC_scMultiome/atac", load = T)

# re-clustering
ao$celltype <- mo$celltype %>% unname()
ao$dpt <- mo$dpt %>% unname()

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
ao <- addImputeWeights(ao)

ao <- addGroupCoverages(
  ArchRProj = ao,
  groupBy = "celltype",
  force = T
)
getGroupBW(
  ArchRProj = ao,
  tileSize = 50,
  groupBy = "celltype"
)
pathToMacs2 <- findMacs2()
ao <- addReproduciblePeakSet(
  ArchRProj = ao,
  groupBy = "celltype",
  pathToMacs2 = pathToMacs2
)

ao <- addPeakMatrix(ao)
ao <- addMotifAnnotations(ao, motifSet = "cisbp", name = "Motif", force = T)
ao <- addBgdPeaks(ao, force = T)
ao <- addDeviationsMatrix(
  ArchRProj = ao,
  peakAnnotation = "cis",
  force = T
)

ao <- addPeak2GeneLinks(
  ArchRProj = ao,
  useMatrix = "GeneExpressionMatrix",
  maxDist = 25e4,
  addEmpiricalPval = T
)
p2g <- getPeak2GeneLinks(
  ArchRProj = ao,
  FDRCutOff = 0.05,
  returnLoops = F
)
p2g <- p2g[p2g$EmpPval < 0.05, ]
saveRDS(p2g, "middata/hADSC_scMultiome/p2g_250k.rds")

# TF 2 gene ---------------------------------------------------------------

p2g <- readRDS("middata/hADSC_scMultiome/p2g_250k.rds")
p2gMeta <- metadata(p2g)

geneNum <- p2gMeta$geneSet$name[p2g$idxRNA] %>% table()
nData <- geneNum %>% as.data.table() %>% set_colnames(c("gene", "N"))
zGene <- p2gMeta$geneSet$name %>% setdiff(nData$gene)
nData <- data.table(gene = zGene, N = 0) %>% rbind(nData)

nData <- nData[order(N)]
nData[, order := 1:.N][]
table(nData$N)

ggplot(nData, aes(x = order, y = N)) +
  geom_point() +
  theme(
    aspect.ratio = 1
  )

nData[, sy := rescale(N)][]
nData[, sx := rescale(order)][]

nData[which.min(1*nData$sy - nData$sx), N]

dorcGene <- nData[N >= 5, gene]
saveRDS(dorcGene, "middata/hADSC_scMultiome/dorcGene.rds")

# peak 2 TF enrichment
bdgPeak <- readRDS("middata/hADSC_scMultiome/atac/Background-Peaks.rds")
bdgPeakMeta <- rowRanges(bdgPeak) %>% as.data.table()

bdgIdx <- assay(bdgPeak)
saveRDS(bdgIdx, "middata/hADSC_scMultiome/bdgIdx.rds")

peakAnno <- getPeakAnnotation(ao, "cis")
motifMat <- readRDS(peakAnno$Matches)
motifPeakMeta <- rowRanges(motifMat) %>% as.data.table()
motifPeakMeta[, id := str_c(seqnames, "_", start, "_", end)]
motifData <- assay(motifMat)

peakMeta <- getPeakSet(ao) %>% as.data.table()
peakMeta[, id := str_c(seqnames, "_", start, "_", end)]
rownames(motifData) <- peakMeta$id

colTF <- colnames(motifData)
motifData <- summary(motifData) %>% as.data.table()
motifData[, j := colTF[j]][]
motifData <- motifData[, .(x = any(x)), by = c("i", "j")]

saveRDS(motifData, "middata/hADSC_scMultiome/ao_motifData.rds")
motifData <- readRDS("middata/hADSC_scMultiome/ao_motifData.rds")

library(progress)
pb <- progress_bar$new(
  format = "running [:bar]:percent | eta :eta | already :elapsed | i = :current",
  total = length(keepGene),
  clear = F,
  stream = stdout())

geneReg <- env()
for(i in keepGene) {
  pb$tick()
  
  corPeak <- p2gData[gene == i, idxATAC]
  corPeakBdg <- bdgIdx[corPeak, ] %>% as.vector() %>% unique()
  
  TFpeak <- motifData[i %in% corPeak, .N, by = j]
  TFbdg <- motifData[i %in% corPeakBdg, .N, by = j]
  
  TFpeak[, Nt := length(corPeak)]
  TFpeak[, B := TFbdg$N[match(TFpeak$j, TFbdg$j)]]
  TFpeak[is.na(B), B := 0]
  TFpeak[, Bt := length(corPeakBdg)]
  
  TFpeak[, p := fisher.test(
    matrix(c(N, Nt, B, Bt), nrow = 2),
    alternative = "greater")$p.value, by = j]
  
  geneReg[[i]] <- TFpeak
}
saveRDS(geneReg, "middata/hADSC_scMultiome/ao_geneReg_full.rds")

geneRegData <- as.list(geneReg)
map_int(geneRegData, nrow) %>% table()
geneRegData <- imap(geneRegData, ~ .x[, gene := .y])

geneTF <- do.call(rbind, geneRegData)
colnames(geneTF)[1] <- "TF"
geneTF <- geneTF[p < 0.05]
geneTF[, TFname := str_remove(TF, ":.*")]

saveRDS(geneTF, "middata/hADSC_scMultiome/ao_geneTF.rds")

# enhancer ----------------------------------------------------------------

# after 'merge all peak' part of CUT&TAG_hADSC_PostAnalysis.R

ao <- addCoAccessibility(
  ArchRProj = ao,
  maxDist = 25e4
)

cA <- getCoAccessibility(
  ao,
  corCutOff = 0.25
)
cA <- cA$CoAccessibility %>% as.data.table()
saveRDS(cA, "ao_cA.rds")

allPeakMeta <- readRDS("allPeakMeta.rds")
genePro <- fread("promoter2k_sort.bed")
genePro %<>% set_colnames(c("chr", "start", "end", "gene"))

e2gList <- list()

pb <- progress_bar$new(
  format = "running [:bar]:percent | eta :eta | already :elapsed | i = :current",
  total = nrow(genePro),
  clear = F,
  stream = stdout())

for(i in 1:nrow(genePro)) {
  pb$tick()
  
  p2k <- c(genePro[i, start], genePro[i, end])
  loop <- cA[seqnames == genePro[i, chr] & (start %b% p2k | end %b% p2k)]
  
  if(nrow(loop) == 0) {
    e2gList[[genePro[i, gene]]] <- "none"
    next
  }
  
  loop$start_side <- "none"
  loop$end_side <- "none"
  loop$side <- -1
  
  loop[start %b% c(genePro[i, start], genePro[i, end]), start_side := "pro"]
  loop[end %b% c(genePro[i, start], genePro[i, end]), end_side := "pro"]
  
  loop[start_side == "none", side := start]
  loop[end_side == "none", side := end]
  
  enhLoop <- loop[side != -1, .(cor = mean(correlation)), by = side]
  
  if(nrow(enhLoop) == 0) {
    e2gList[[genePro[i, gene]]] <- "none"
    next
  }
  
  for(j in 1:nrow(enhLoop)) {
    
    pos <- enhLoop[j, side]
    overPeak <- allPeakMeta[chr == genePro[i, chr] & start <= pos & end >= pos]
    overPeak <- overPeak[H3K27ac != "-" | H3K4me1 != "-"]
    overPeak <- overPeak[!(str_detect(gene, genePro[i, gene]) & type == "pro")]
    
    if(nrow(overPeak) == 0) {enhLoop[j, peak := "no"]} else {
      enhLoop[j, peak := overPeak$id]
    }
  }
  
  enhPeak <- enhLoop[peak != "no", .(cor = mean(cor)), by = peak]
  
  if(nrow(enhPeak) == 0) {
    e2gList[[genePro[i, gene]]] <- "none"
    next
  }
  
  e2gList[[genePro[i, gene]]] <- enhPeak
}

saveRDS(e2gList, "e2gList.rds")
