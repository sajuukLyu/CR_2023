
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu, Cheng Li Lab, Peking University
# email: lvyulin@pku.edu.cn
# date: 2023.4.24
#
# ---

# load package ------------------------------------------------------------

library(tidybayes)
library(data.table)
library(glue)

source("R/utils.R")

# after peak calling process ----------------------------------------------

# merge peak
fileName <- list.files("peak", "broadPeak")
sampName <- str_remove(fileName, "_peak.*")

usedHis <- c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac")

sampList <- map(usedHis, ~ str_subset(sampName, .x)) %>% set_names(usedHis)

allPeakList <- list()
for(i in names(sampList)) {
  message(i)
  allPeakList[[i]] <- sampList[[i]] %>% str_c("peak/", ., "_peaks.broadPeak") %>% map(~ fread(.x)) %>% do.call(rbind, .)
}

iwalk(allPeakList, ~ {fwrite(.x[, 1:3], str_c("merge/", .y, "_all.bed"), sep = "\t", col.names = F)})

mergeDir <- "merge"
mergeCMD <- glue(
  "cat {mergeDir}/{names(sampList)}_all.bed | sort -k1,1 -k2,2n > {mergeDir}/{names(sampList)}_all.sort.bed
  bedtools merge -d 0 -i {mergeDir}/{names(sampList)}_all.sort.bed > {mergeDir}/{names(sampList)}_merge.bed")
mergeCMD
write.table(c("#!/bin/bash", mergeCMD), glue("{codeDir}/{mergeDir}.sh"), row.names = F, col.names = F, quote = F, sep = "\n")

overCMD <- c()
for(i in names(sampList)) {
  overCMD <- glue(
    "bedtools intersect -wa -wb -a {mergeDir}/{i}_merge.bed -b peak/{sampList[[i]]}_peaks.broadPeak > \\
    {mergeDir}/{sampList[[i]]}_over.bed") %>% c(overCMD, .)
}
overCMD[1]
dir.create(glue("{codeDir}/{mergeDir}"))
setCMD(overCMD, glue("{codeDir}/{mergeDir}"), 6)

allPeakData <- map(names(sampList), ~ fread(str_c(mergeDir, "/", .x, "_merge.bed")) %>% set_colnames(c("chr", "start", "end"))) %>% set_names(names(sampList))
allPeakData <- map(allPeakData, ~ {.x[, id := str_c(chr, "_", start, "_", end)]})
saveRDS(allPeakData, "allPeakData.rds")
allPeakData <- readRDS("allPeakData.rds")

hasPeakList <- list()
for(i in names(sampList)) {
  message(i)
  overData <- map(sampList[[i]], ~ fread(str_c(mergeDir, "/", .x, "_over.bed"))) %>% set_names(sampList[[i]])
  overData <- map(overData, ~ {.x[, id := str_c(V1, "_", V2, "_", V3)]})
  
  hasPeakList[[i]] <- map(overData, ~ {allPeakData[[i]]$id %in% .x$id}) %>% do.call(cbind, .) %>% set_rownames(allPeakData[[i]]$id)
}
saveRDS(hasPeakList, "hasPeakList.rds")
hasPeakList <- readRDS("hasPeakList.rds")


# z score
fileName <- list.files("spbw", ".bw") # .bw files for each repeat

usedHis <- c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac")

bwList <- map(usedHis, ~ str_subset(fileName, .x)) %>% set_names(usedHis)

binCMD <- c()
for(i in names(bwList)) {
  message(i)
  binCMD <- glue(
    "{deepDir}/multiBigwigSummary bins -b {bws} -bs 1000 -p 24 \\
    -o bin/{i}.npz --outRawCounts bin/{i}.tsv",
    bws = str_c("spbw/", bwList[[i]]) %>% str_c(collapse = " ")) %>% c(binCMD, .)
}
cat(binCMD[1])

dir.create("bin")
dir.create(glue("{codeDir}/bin"))
setCMD(binCMD, glue("{codeDir}/bin"), 6)

meanList <- list()
sdList <- list()

for(i in names(bwList)) {
  message(i)
  binData <- fread(str_c("bin/", i, ".tsv"))
  binData <- binData[, -1:-3] %>% as.matrix()
  binData[is.nan(binData)] <- 0
  meanList[[i]] <- colMeans(binData)
  sdList[[i]] <- apply(binData, 2, sd)
}

saveRDS(meanList, "meanList.rds")
saveRDS(sdList, "sdList.rds")

for(i in names(bwList)) {
  message(i)
  
  bdgFile <- bwList[[i]] %>% str_replace(".bw", "_treat_pileup.bdg")
  outFile <- bwList[[i]] %>% str_replace(".bw", "_zscore.bdg")
  
  bdgData <- map(bdgFile, ~ fread(str_c("bdg/", .x)))
  
  bdgData <- map2(bdgData, meanList[[i]], ~ .x[, V4 := V4 - .y])
  bdgData <- map2(bdgData, sdList[[i]], ~ .x[, V4 := V4 / .y])
  
  walk2(bdgData, outFile, ~ fwrite(.x, str_c("zbdg/", .y), sep = "\t", col.names = F))
}

dir.create(glue("{codeDir}/zbw"))

zbdgFile <- list.files("zbdg")
zbwFile <- str_replace(zbdgFile, ".bdg", ".bw")

bdg2bwDir = "bedGraphToBigWig"
refLen = "genome.fa.fai"

zbwCMD <- glue(
  "{bdg2bwDir} zbdg/{zbdgFile} {refLen} zbw/{zbwFile}")
cat(zbwCMD[1])

setCMD(zbwCMD, glue("{codeDir}/zbw"), 6)


# group z score
fileName <- list.files("bw", ".bw")

usedHis <- c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac")

bwList <- map(usedHis, ~ str_subset(fileName, .x)) %>% set_names(usedHis)

binCMD <- c()
for(i in names(bwList)) {
  message(i)
  binCMD <- glue(
    "{deepDir}/multiBigwigSummary bins -b {bws} -bs 1000 -p 24 \\
    -o mbin/{i}.npz --outRawCounts mbin/{i}.tsv",
    bws = str_c("bw/", bwList[[i]]) %>% str_c(collapse = " ")) %>% c(binCMD, .)
}
cat(binCMD[1])

dir.create("mbin")
dir.create(glue("{codeDir}/mbin"))
setCMD(binCMD, glue("{codeDir}/mbin"), 6)

meanList <- list()
sdList <- list()

for(i in names(bwList)) {
  message(i)
  binData <- fread(str_c("mbin/", i, ".tsv"))
  binData <- binData[, -1:-3] %>% as.matrix()
  binData[is.nan(binData)] <- 0
  meanList[[i]] <- colMeans(binData)
  sdList[[i]] <- apply(binData, 2, sd)
}

saveRDS(meanList, "meanListMerge.rds")
saveRDS(sdList, "sdListMerge.rds")

for(i in names(bwList)) {
  message(i)
  
  bdgFile <- bwList[[i]] %>% str_replace(".bw", "_treat_pileup.bdg")
  outFile <- bwList[[i]] %>% str_replace(".bw", "_zscore.bdg")
  
  bdgData <- map(bdgFile, ~ fread(str_c("mbdg/", .x)))
  
  bdgData <- map2(bdgData, meanList[[i]], ~ .x[, V4 := V4 - .y])
  bdgData <- map2(bdgData, sdList[[i]], ~ .x[, V4 := V4 / .y])
  
  walk2(bdgData, outFile, ~ fwrite(.x, str_c("mzbdg/", .y), sep = "\t", col.names = F))
}

dir.create(glue("{codeDir}/mzbw"))

zbdgFile <- list.files("mzbdg")
zbwFile <- str_replace(zbdgFile, ".bdg", ".bw")

zbwCMD <- glue(
  "{bdg2bwDir} mzbdg/{zbdgFile} {refLen} mzbw/{zbwFile}")
cat(zbwCMD[1])

setCMD(zbwCMD, glue("{codeDir}/mzbw"), 6)


# value of merged peak
usedHis <- c("H3K4me3", "H3K4me1", "H3K27me3", "H3K27ac")

fileName <- list.files("zbw", ".bw")

bwList <- map(usedHis, ~ str_subset(fileName, .x)) %>% set_names(usedHis)
saveRDS(bwList, "bwList.rds")
bwList <- readRDS("bwList.rds")

scoreCMD <- c()
for(i in names(bwList)) {
  bwFile <- bwList[[i]]
  bedFile <- glue("merge/{i}_merge.bed")
  scoreCMD <- glue(
    "{deepDir}/computeMatrix scale-regions -S {bws} -R {bedFile} \\
    -m 1000 -a 2000 -b 2000 -bs 20 -p 24 \\
    -o score/{i}.mat --outFileNameMatrix score/{i}.tsv --outFileSortedRegions score/{i}.bed \\\n
    {deepDir}/plotHeatmap -m score/{i}.mat --colorMap RdBu_r \\
    -o score/{i}.heatmap.pdf --outFileSortedRegions score/{i}.sort.bed",
    bws = str_c("zbw/", bwFile) %>% str_c(collapse = " ")) %>% c(scoreCMD, .)
}
cat(scoreCMD[1])
dir.create("score")
dir.create("code/score")
setCMD(scoreCMD, glue("{codeDir}/score"), 6)

allPeakData <- readRDS("allPeakData.rds")

valList <- list()
for(i in names(bwList)) {
  message(i)
  allVal <- fread(str_c("score/", i, ".tsv"), skip = 3, header = F)
  peakMeta <- fread(str_c("score/", i, ".bed"))
  
  colnames(peakMeta)[1] <- "chr"
  peakMeta[, id := str_c(chr, "_", start, "_", end)]
  
  allVal <- allVal[match(allPeakData[[i]]$id, peakMeta$id)]
  allVal <- as.matrix(allVal)
  allVal[is.na(allVal)] <- 0
  
  valMat <- allVal %>% apply(1, function(x) tapply(x, rep(bwList[[i]], each = 250), function(y) mean(y[101:150]))) %>% t
  valMat %<>% set_rownames(allPeakData[[i]]$id)
  
  valList[[i]] <- valMat
}

saveRDS(valList, "valList.rds")

