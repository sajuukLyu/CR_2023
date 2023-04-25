
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

source("R/utils.R")

`%ni%` <- function(a, b) return(! a %in% b)

# peak annotate -----------------------------------------------------------

setwd("middata/hADSC_CUT&TAG")

allPeakData <- readRDS("allPeakData.rds")
valList <- readRDS("valList.rds")
hasPeakList <- readRDS("hasPeakList.rds")

# gene overlap
allPeakData <- map(allPeakData, ~ .x[, chr :=  str_c("chr", chr)])

iwalk(allPeakData, ~ {fwrite(.x[, 1:4], str_c(.y, "_data.tsv"), sep = "\t", row.names = F, col.names = F)})

geneCMD <- str_c(
  "cat ", names(allPeakData), "_data.tsv | sort -k1,1 -k2,2n > ", names(allPeakData), "_data_sort.tsv\n",
  "bedtools closest -a ", names(allPeakData), "_data_sort.tsv -b promoter2k_sort.bed -d > ", names(allPeakData), "_p2k.bed"
)
cat(geneCMD)
walk(geneCMD, system)

for(i in names(allPeakData)) {
  message(i)
  peakOver <- fread(str_c(i, "_p2k.bed"))
  peakOver <- peakOver[, .(gene = str_c(V8, collapse = "|"), dist = str_c(V9, collapse = "|")), by = V4]
  
  allPeakData[[i]]$gene <- peakOver[match(allPeakData[[i]]$id, V4), gene]
  allPeakData[[i]]$dist <- peakOver[match(allPeakData[[i]]$id, V4), dist]
  
  allPeakData[[i]][, type := "enh"]
  allPeakData[[i]][dist %>% str_split("\\|") %>% map_lgl(~ sum(as.numeric(.x)) == 0), type := "pro"]
}

# filter peaks
res <- map2(hasPeakList, valList, ~ roc(as.vector(.x[, rep(1:ncol(.x), each = 2)]), as.vector(.y)))
th <- map(res, ~ .x$thresholds[which.max(.x$sensitivities + .x$specificities - 1)])

saveRDS(th, "th.rds")
th <- readRDS("th.rds")

for(i in names(valList)) {
  message(i)
  passTh <- valList[[i]] %>% apply(1, function(x) max(x)) > th[[i]]
  allPeakData[[i]]$passTh <- passTh
}

scaleList <- map(valList, ~ {apply(.x, 1, function(x) {x/max(x)}) %>% t})
repList <- map(scaleList, ~ {apply(.x, 1, function(x) tapply(x, str_remove(colnames(.x), "_H3.*"), function(y) abs(y[1] - y[2]) < 0.5)) %>% t})

for(i in names(repList)) {
  message(i)
  passRep <- apply(repList[[i]], 1, all)
  allPeakData[[i]]$passRep <- passRep
}

peakData <- map(allPeakData, ~ .x[chr %ni% c("chrX", "chrY") & passTh & passRep])
saveRDS(peakData, "peakData.rds")

# global change -----------------------------------------------------------

usedList <- map2(valList, peakData, ~ .x[.y$id, ])
avgList <- map(usedList, ~ {apply(.x, 1, function(x) tapply(x, str_remove(colnames(.x), "_H.*") %>% {factor(., unique(.))}, mean)) %>% t})

diffData <- avgList$H3K4me3 %>% {. - .[, 1]}
diffData <- avgList$H3K27me3 %>% {. - .[, 1]}
diffData <- avgList$H3K4me1 %>% {. - .[, 1]}
diffData <- avgList$H3K27ac %>% {. - .[, 1]}
colSums(diffData)

plotData <- diffData %>% as.data.table(T) %>% melt(id.vars = "rn", variable.name = "samp", value.name = "val")
plotData$samp %<>% factor(colnames(diffData))

ggplot(plotData, aes(x = samp, y = val)) +
  geom_violin(scale = "width") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plotData[, type := "un"]
plotData[val > 1, type := "up_1~2"]
plotData[val > 2, type := "up_2~5"]
plotData[val > 5, type := "up_>5"]
plotData[val < -1, type := "down_1~2"]
plotData[val < -2, type := "down_2~5"]
plotData[val < -5, type := "down_>5"]

pctData <- plotData[, .(n = .N), by = c("samp", "type")]
pctData <- pctData[type != "un"]
pctData[str_detect(type, "down"), n := -n][]
pctData[, n := n/nrow(diffData)]
pctData$samp %<>% factor(colnames(diffData))

changeCol <- structure(
  ArchR::ArchRPalettes$solarExtra[c(9:7, 1:3)] %>% unname(),
  names = c("up_>5", "up_2~5", "up_1~2", "down_>5", "down_2~5", "down_1~2")
)
pctData$type %<>% factor(names(changeCol))

sumData <- pctData[, .(diff = sum(n[str_detect(type, "up")]) + sum(n[str_detect(type, "down")])), by = samp]

ggplot(pctData, aes(x = samp, y = n * 100)) +
  geom_col(aes(fill = type), position = "stack") +
  geom_hline(yintercept = 0, size = .5, color = "gray50") +
  geom_point(data = sumData, aes(y = diff * 100)) +
  geom_path(data = sumData[1:6], group = "all", aes(y = diff * 100)) +
  scale_fill_manual(values = changeCol) +
  labs(x = "", y = str_c("pct. of merged peaks\n(n = ", nrow(diffData), ")"), fill = "z-score\nchange") +
  theme(
    aspect.ratio = 2/3,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# cluster peak ------------------------------------------------------------

i = "H3K4me3"

usedList <- map2(scaleList, peakData, ~ .x[.y$id, ])

heatData <- usedList[[i]]

colorFun <- colorRamp2(seq(1, 0, len = 9), brewer.pal(9, "RdBu"))
p <- Heatmap(
  heatData, col = colorFun, border = F,
  cluster_rows = T, cluster_columns = F,
  show_row_names = F, show_column_names = T,
  show_row_dend = F, show_column_dend = F,
  row_km = 24, row_km_repeats = 100)

peakOrder <- row_order(p)
peakType <- data.table(
  type = rep(names(peakOrder), map_int(peakOrder, length)),
  peak = unlist(peakOrder) %>% unname()
)

png(glue("{i}_clus.png"), width = 10, height = 14, units = "in", res = 300)
Heatmap(
  heatData[peakType$peak, ], col = colorFun, border = T,
  cluster_rows = F, cluster_columns = F,
  show_row_names = F, show_column_names = T,
  show_row_dend = F, show_column_dend = F,
  column_names_rot = 45,
  row_split = peakType$type,
  heatmap_legend_param = list(
    title = "Scaled expression",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(6.5, "in"),
  height = unit(10.5, "in"),
  use_raster = T
)
dev.off()

peakData[[i]][peakType$peak, cluster := peakType$type][]
peakData[[i]]$group <- "un"
peakData[[i]][cluster %in% c(), group := "a1"]
peakData[[i]][cluster %in% c(), group := "a2"]
peakData[[i]][cluster %in% c(), group := "a3"]
peakData[[i]][cluster %in% c(), group := "a4"]
peakData[[i]][cluster %in% c(), group := "a5"]
table(peakData[[i]]$group)

png(glue("{i}.heatmap.png"), width = 7, height = 8, units = "in", res = 300)
Heatmap(
  heatData[peakData[[i]]$id, ], col = colorFun, border = F,
  cluster_rows = F, cluster_columns = F,
  show_row_names = F, show_column_names = T,
  show_row_dend = F, show_column_dend = F,
  column_names_rot = 45,
  row_split = peakData$H3K4me3$group,
  heatmap_legend_param = list(
    title = "Scaled expression",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(4, "in"),
  height = unit(5.5, "in"))
dev.off()

geneList <- tapply(peakData[[i]][type == "pro", gene], peakData[[i]][type == "pro", group], c) %>%
  map(~ {str_c(.x, collapse = "|") %>% str_split("\\|") %>% pluck(1) %>% unique})
geneList$un <- NULL

saveRDS(geneList, glue("{i}.geneList.rds"))

geneMeta <- data.table(
  gene = geneList %>% unlist,
  cluster = rep(names(geneList), map_int(geneList, length))
)
fwrite(geneMeta, glue("{i}_clus_geneList.csv"), sep = ",")


bwList <- readRDS("bwList.rds")
scoreData <- fread(glue("score/{i}.tsv"), skip = 3, header = F)

peakMeta <- fread(glue("score/{i}.bed"))
colnames(peakMeta)[1] <- "chr"
peakMeta[, id := str_c(chr, "_", start, "_", end)][]

dim(scoreData)
scoreData <- scoreData[match(peakData[[i]]$id, peakMeta$id)]
scoreData <- as.matrix(scoreData)
scoreData[is.na(scoreData)] <- 0

posOrder <- map2(
  (1:length(bwList[[i]]) * 250 - 249)[match(colnames(heatData), bwList[[i]])],
  (1:length(bwList[[i]]) * 250)[match(colnames(heatData), bwList[[i]])], ~ .x:.y) %>% do.call(c, .)
scoreData <- scoreData[, posOrder]
peakData[[i]]$rs <- rowSums(scoreData)
peakOrder <- peakData[[i]][group != "un"][order(group, -rs)]
scoreData <- scoreData[match(peakOrder$id, peakData[[i]]$id), ]
sampOrder <- colnames(heatData) %>% str_remove_all(".bw") %>% rep(each = 250) %>% {factor(., unique(.))}

valData <- valList[[i]][peakOrder$id, ]
rmax <- apply(valData, 1, max)
scoreData <- sweep(scoreData, 1, rmax, "/")

colorFun <- colorRamp2(seq(0, 1, len = 9)^0.5, brewer.pal(9, "Reds"))

markCol <- c(
  getMatColor("red", 2)[2],
  getMatColor("deep-orange", 1),
  getMatColor("orange", 1),
  getMatColor("amber", 1),
  getMatColor("yellow", 1))
markCol %>% show_col()

pdf(glue("{i}.heatmap.pdf"), width = 10, height = 12)
Heatmap(
  scoreData, col = colorFun, border = F,
  cluster_rows = F, cluster_columns = F,
  show_row_names = F, show_column_names = F,
  show_row_dend = F, show_column_dend = F,
  row_split = peakOrder$group,
  row_gap = unit(0.05, "in"),
  column_split = sampOrder,
  column_title_rot = 90,
  column_gap = unit(c(rep(c(0, 0.05), length(bwList[[i]])/2-1), 0), "in"),
  right_annotation = rowAnnotation(
    type = anno_block(gp = gpar(fill = markCol, col = NA), width = unit(0.2, "in"))
  ),
  heatmap_legend_param = list(
    title = "Relative Normalized Signal",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(8, "in"),
  height = unit(8, "in"))
dev.off()

# dynamic pca -------------------------------------------------------------

dynData <- map2(valList, peakData, ~ {.x[.y[group %ni% c("un"), id], ]})

PCAdata <- map(dynData, ~ prcomp_irlba(.x, n = 3, scale. = F))
PCprop <- map(PCAdata, ~ {(.x$sdev)^2 %>% {round(. / sum(.), 3) * 100}})

plotData <- imap(PCAdata, ~ as.data.table(
  .x$rotation)[, samp := str_remove_all(colnames(dynData$H3K4me3), "_H3.*")][, mark := .y]) %>% do.call(rbind, .)

labData <- plotData[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = c("samp", "mark")]

sampColor <- structure(c(
  getMatColor("blue", 1),
  getMatColor("green", 2),
  getMatColor("yellow", 1),
  getMatColor("amber", 1),
  getMatColor("orange", 1),
  getMatColor("deep-orange", 1),
  getMatColor("red", 2),
  getMatColor("purple", 2)), names = c(
    "ADSC", "S1D4", "S1D8",
    "S2D4", "S2D8", "S2D12", "S2D16",
    "IPSno3", "IPSno7", "H1", "H9"))
sampColor %>% show_col()

i = "H3K4me3"
i = "H3K27me3"
i = "H3K4me1"
i = "H3K27ac"

ggplot(plotData[mark == i], aes(x = PC1, y = PC2)) +
  geom_point(aes(color = samp), size = 3, show.legend = T) +
  geom_text_repel(data = labData[mark == i], aes(label = samp), size = 2) +
  scale_color_manual(values = sampColor) +
  labs(
    x = str_c("PC1 (", PCprop[[i]][1], "%)"),
    y = str_c("PC2 (", PCprop[[i]][2], "%)"),
    title = i, color = "") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_blank()
  )

ggsave(str_c(i, ".PCA.pdf"), width = 4, height = 3)

# dynamic corr ------------------------------------------------------------

dynData <- map2(valList, peakData, ~ {.x[.y[group %ni% c("un"), id], ]})

avgData <- map(dynData, ~ {apply(.x, 1, function(x) tapply(x, str_remove_all(colnames(dynData$H3K4me3), "_H3.*"), mean)) %>% t})
sampOrder <- str_remove_all(colnames(dynData$H3K4me3), "_H3.*") %>% unique()
corData <- map(avgData, ~ cor(.x[, sampOrder], method = "pearson"))

disp <- seq(0, 1, len = 9)

colorFun <- colorRamp2(disp, ArchR::ArchRPalettes$solarExtra)

i <- "H3K4me3"
i <- "H3K27me3"
i <- "H3K4me1"
i <- "H3K27ac"

pdf(str_c(i, ".cor.pdf"), width = 5, height = 5)
Heatmap(
  corData[[i]], col = colorFun,
  cluster_rows = F, cluster_columns = F,
  column_title = i,
  heatmap_legend_param = list(
    title = "Pearson correlation",
    title_position = "lefttop-rot",
    legend_height = unit(1.5, "in")),
  width = unit(3, "in"),
  height = unit(3, "in"))
dev.off()

# merge all peak ----------------------------------------------------------

map(peakData, ~ .x[passTh & passRep][, 1:3]) %>% do.call(rbind, .) %>%
  fwrite("allPeak.tsv", sep = "\t", row.names = F, col.names = F)

"cat allPeak.tsv | sort -k1,1 -k2,2n > allPeak_sort.tsv"
"bedtools merge -d 0 -i allPeak_sort.tsv > allPeak_merge.tsv"

allPeakMeta <- fread("allPeak_merge.tsv")
allPeakMeta[, id := str_c("peak_", 1:.N)]
fwrite(allPeakMeta, "allPeak_merge.tsv", sep = "\t", row.names = F, col.names = F)
colnames(allPeakMeta) <- c("chr", "start", "end", "id")

allPeakMeta[, start := start - 48e3]
allPeakMeta[, end := end + 48e3]
allPeakMeta[start < 0, start := 0]
fwrite(allPeakMeta, "allPeak_merge_48k.tsv", sep = "\t", row.names = F, col.names = F)

markPeak <- c("H3K4me3_data_sort.tsv", "H3K4me1_data_sort.tsv", "H3K27ac_data_sort.tsv", "H3K27me3_data_sort.tsv")
overCMD <- str_c(
  "bedtools intersect -wb -a allPeak_merge.tsv",
  " -b ", markPeak, " > ", str_replace(markPeak, "_sort.tsv", "_overPeak.tsv"))
write.table(c("#!/bin/bash", overCMD), "mark_over.sh", sep = "\n", row.names = F, col.names = F, quote = F)

for(i in names(peakData)) {
  message(i)
  
  peakOver <- fread(str_c(i, "_data_overPeak.tsv"))
  peakOver <- peakOver[, .(peak = str_c(V8, collapse = ";")), by = V4]
  allPeakMeta[[i]] <- "-"
  
  allPeakMeta[[i]][match(peakOver$V4, allPeakMeta$id)] <- peakOver$peak
}

"bedtools closest -a allPeak_merge.tsv -b promoter2k_sort.bed -d > allPeak_p2k.bed"
"bedtools intersect -wb -a allPeak_merge_48k.tsv -b /mnt/d/project/1_big_scRNAseq/P1/promoter2k_sort.bed > allPeak_p50k.bed"

peakOver <- fread("allPeak_p2k.bed")
peakOver <- peakOver[, .(gene = str_c(V8, collapse = "|"), dist = str_c(V9, collapse = "|")), by = V4]

allPeakMeta$gene <- peakOver[match(allPeakMeta$id, V4), gene]
allPeakMeta$dist <- peakOver[match(allPeakMeta$id, V4), dist]

allPeakMeta$type <- "enh"
allPeakMeta[dist %>% str_split("\\|") %>% map_lgl(~ sum(as.numeric(.x)) == 0), type := "pro"]
table(allPeakMeta$type)
allPeakMeta[, l := end - start]
allPeakMeta$l %>% summary()

names(valList)
overThList <- map2(valList, th, ~ .x > .y)
cn <- overThList$H3K4me3 %>% colnames() %>% str_remove_all("_H.*")
overThList <- map(overThList, ~ {apply(.x, 1, function(x) tapply(x, cn, all)) %>% t})
overThList <- map(overThList, ~ .x[, unique(cn)])
saveRDS(overThList, "overThList.rds")

sigMatList <- list()
for(i in names(overThList)) {
  message(i)
  peakSp <- allPeakMeta[[i]] %>% str_split(";") %>% map(~ setdiff(.x, "-"))
  peakVal <- map(peakSp, ~ {
    if(length(.x) == 0) {rep(0, ncol(overThList[[i]]))} else
      if(length(.x) == 1) {overThList[[i]][.x, ]} else {
        apply(overThList[[i]][.x, ], 2, any)
      }
  })
  
  sigMatList[[i]] <- do.call(rbind, peakVal) %>% set_rownames(allPeakMeta$id)
}
saveRDS(sigMatList, "sigMatList_overTh.rds")
sigMatList <- readRDS("sigMatList_overTh.rds")

sampOrder <- colnames(sigMatList$H3K4me3) %>% str_remove_all("_H3.*")

sigMatList <- map(sigMatList, ~ {set_colnames(.x, sampOrder) %>% as.data.table(T)})

for(i in sampOrder) {
  message(i)
  allPeakMeta[[str_c(i, "_pro")]] <- "NoSig"
  
  allPeakMeta[[str_c(i, "_pro")]][sigMatList$H3K4me3[[i]] & !sigMatList$H3K27me3[[i]]] <- "Act"
  allPeakMeta[[str_c(i, "_pro")]][sigMatList$H3K4me3[[i]] & sigMatList$H3K27me3[[i]]] <- "Biv"
  allPeakMeta[[str_c(i, "_pro")]][!sigMatList$H3K4me3[[i]] & sigMatList$H3K27me3[[i]]] <- "Ina"
}

for(i in sampOrder) {
  message(i)
  allPeakMeta[[str_c(i, "_enh")]] <- "NoSig"
  
  allPeakMeta[[str_c(i, "_enh")]][sigMatList$H3K4me1[[i]] & !sigMatList$H3K27ac[[i]]] <- "Poi"
  allPeakMeta[[str_c(i, "_enh")]][sigMatList$H3K4me1[[i]] & sigMatList$H3K27ac[[i]]] <- "Act"
  allPeakMeta[[str_c(i, "_enh")]][!sigMatList$H3K4me1[[i]] & sigMatList$H3K27ac[[i]]] <- "Nor"
}

usedC <- str_subset(colnames(allPeakMeta), "pro")
proType <- allPeakMeta[, lapply(.SD, table), .SDcols = usedC]
usedC <- str_subset(colnames(allPeakMeta), "enh")
enhType <- allPeakMeta[, lapply(.SD, table), .SDcols = usedC]

for(i in names(peakData)) {
  message(i)
  peakSp <- allPeakMeta[[i]] %>% str_split(";") %>% map(~ setdiff(.x, "-"))
  peakType <- map_chr(peakSp, ~ {
    if(length(.x) == 0) {"-"} else
      if (.x %ni% peakData[[i]]$id) {"-"} else{
        peakData[[i]][match(.x, id), group] %>% unique %>% sort %>% str_c(collapse = "|")
      }})
  allPeakMeta[[str_c(i, "_group")]] <- peakType
}

allPeakMeta$H3K4me3_group %>% table
allPeakMeta[H3K4me3_group %>% str_detect("\\|un"), H3K4me3_group := str_remove(H3K4me3_group, ".un")]
# allPeakMeta[H3K4me3_group %>% str_detect("a1\\|"), H3K4me3_group := "a1"]
# allPeakMeta[H3K4me3_group == "a2|a3", H3K4me3_group := "a3"]
# allPeakMeta[H3K4me3_group == "a2|a4", H3K4me3_group := "a4"]
# allPeakMeta[H3K4me3_group == "a2|a5", H3K4me3_group := "a2"]
# allPeakMeta[H3K4me3_group == "a3|a4", H3K4me3_group := "a3"]
# allPeakMeta[H3K4me3_group == "a3|a5", H3K4me3_group := "a3"]
# allPeakMeta[H3K4me3_group == "a4|a5", H3K4me3_group := "a4"]

allPeakMeta$H3K27me3_group %>% table
allPeakMeta[H3K27me3_group %>% str_detect("\\|un"), H3K27me3_group := str_remove(H3K27me3_group, ".un")]
# allPeakMeta[H3K27me3_group == "b1|b2", H3K27me3_group := "b1"]
# allPeakMeta[H3K27me3_group == "b1|b3", H3K27me3_group := "b1"]
# allPeakMeta[H3K27me3_group == "b2|b3", H3K27me3_group := "b2"]

allPeakMeta$H3K4me1_group %>% table
allPeakMeta[H3K4me1_group %>% str_detect("\\|un"), H3K4me1_group := str_remove(H3K4me1_group, ".un")]
# allPeakMeta[H3K4me1_group %>% str_detect("c1\\|"), H3K4me1_group := str_remove(H3K4me1_group, "c1.")]
# allPeakMeta[H3K4me1_group == "c2|c3", H3K4me1_group := "c3"]
# allPeakMeta[H3K4me1_group == "c2|c4", H3K4me1_group := "c4"]
# allPeakMeta[H3K4me1_group == "c2|c5", H3K4me1_group := "c2"]
# allPeakMeta[H3K4me1_group == "c3|c4", H3K4me1_group := "c4"]
# allPeakMeta[H3K4me1_group == "c3|c5", H3K4me1_group := "c3"]
# allPeakMeta[H3K4me1_group == "c4|c5", H3K4me1_group := "c4"]

allPeakMeta$H3K27ac_group %>% table
allPeakMeta[H3K27ac_group %>% str_detect("\\|un"), H3K27ac_group := str_remove(H3K27ac_group, ".un")]
# allPeakMeta[H3K27ac_group %>% str_detect("d1\\|"), H3K27ac_group := str_remove(H3K27ac_group, "d1.")]
# allPeakMeta[H3K27ac_group %>% str_detect("\\|d5"), H3K27ac_group := str_remove(H3K27ac_group, ".d5")]
# allPeakMeta[H3K27ac_group == "d2|d3", H3K27ac_group := "d3"]
# allPeakMeta[H3K27ac_group == "d2|d3|d4", H3K27ac_group := "d3"]
# allPeakMeta[H3K27ac_group == "d2|d4", H3K27ac_group := "d4"]
# allPeakMeta[H3K27ac_group == "d3|d4", H3K27ac_group := "d3"]

table(allPeakMeta[, c("H3K4me3_group", "H3K27me3_group")])
table(allPeakMeta[, c("H3K4me1_group", "H3K27ac_group")])

saveRDS(allPeakMeta, "allPeakMeta.rds")

# diff enhancer -----------------------------------------------------------

nData <- allPeakMeta[type == "enh" & H3K4me1_group != "-" & H3K27ac_group != "-", .N, by = c("ADSC_enh", "S2D16_enh")]
nData$ADSC_enh %<>% factor(c("NoSig", "Nor", "Poi", "Act") %>% rev)
nData$S2D16_enh %<>% factor(c("NoSig", "Nor", "Poi", "Act"))

nData$type <- "un"
nData[ADSC_enh == "Act" & S2D16_enh %in% c("NoSig", "Nor", "Poi"), type := "down"]
nData[ADSC_enh %in% c("Poi", "Nor") & S2D16_enh %in% c("NoSig"), type := "down"]
nData[S2D16_enh == "Act" & ADSC_enh %in% c("NoSig", "Nor", "Poi"), type := "up"]
nData[S2D16_enh %in% c("Poi", "Nor") & ADSC_enh %in% c("NoSig"), type := "up"]

S2col <- structure(
  c(getMatColor("red", 1),
    getMatColor("blue", 1),
    "gray80"),
    names = c("up", "down", "un")
)

ggplot(nData, aes(x = S2D16_enh, y = ADSC_enh)) +
  geom_tile(aes(fill = type), color = "black", size = .5, show.legend = F) +
  geom_text(aes(label = N)) +
  scale_fill_manual(values = S2col) +
  coord_fixed() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("enh_N.pdf", width = 3, height = 3)

S2enh <- allPeakMeta[
  type == "enh" & H3K4me1_group != "-" & H3K27ac_group != "-" &
    S2D16_enh == "Act" & ADSC_enh != "Act", id]

ADSCenh <- allPeakMeta[
  type == "enh" & H3K4me1_group != "-" & H3K27ac_group != "-" &
    ADSC_enh == "Act" & S2D16_enh != "Act", id]

geneOver <- fread("allPeak_p50k.bed")
geneOver <- geneOver[V4 %in% c(S2enh, ADSCenh)]
geneOver$type <- "ADSC"
geneOver[V4 %in% S2enh, type := "S2"][]
table(geneOver$type)

geneOver[V8 == "SALL4"]

geneEnh <- geneOver[, .N, by = c("V8", "type")] %>% set_colnames(c("gene", "enh", "N"))
geneData <- dcast(geneEnh, gene ~ enh)
geneData[is.na(geneData)] <- 0
geneData$type <- "un"
geneData[ADSC - S2 <= -4, type := "up"]
geneData[ADSC - S2 >= 4, type := "down"]

pairData <- geneData[, .N, by = c("ADSC", "S2", "type")]

ggplot(pairData, aes(x = S2, y = ADSC)) +
  geom_abline(slope = 1, intercept = 3.5, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -3.5, linetype = "dashed") +
  geom_point(aes(size = log10(N), color = type)) +
  scale_radius() +
  scale_color_manual(values = S2col) +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.key = element_blank()
  )
ggsave("gene_enh_pair_N.pdf", width = 5, height = 3.2)

geneList <- list(
  S2 = geneData[type == "up", gene],
  ADSC = geneData[type == "down", gene]
)

# diff promoter -----------------------------------------------------------

table(allPeakMeta$type)
nData <- allPeakMeta[type == "pro" & (H3K4me3_group != "-" | H3K27me3_group != "-"), .N, by = c("ADSC_pro", "S2D16_pro")]
nData <- rbind(nData, data.table(ADSC_pro = "Biv", S2D16_pro = "NoSig", N = 0))
nData$ADSC_pro %<>% factor(c("Ina", "NoSig", "Biv", "Act") %>% rev)
nData$S2D16_pro %<>% factor(c("Ina", "NoSig", "Biv", "Act"))

nData$type <- "un"
nData[ADSC_pro == "Act" & S2D16_pro %in% c("NoSig", "Ina", "Biv"), type := "down"]
nData[ADSC_pro %in% c("Biv") & S2D16_pro %in% c("Ina", "NoSig"), type := "down"]
nData[ADSC_pro %in% c("NoSig") & S2D16_pro %in% c("Ina"), type := "down"]
nData[S2D16_pro == "Act" & ADSC_pro %in% c("NoSig", "Ina", "Biv"), type := "up"]
nData[S2D16_pro %in% c("Biv") & ADSC_pro %in% c("Ina", "NoSig"), type := "up"]
nData[S2D16_pro %in% c("NoSig") & ADSC_pro %in% c("Ina"), type := "up"]

S2col <- structure(
  c(getMatColor("red", 1),
    getMatColor("blue", 1),
    "gray80"),
    names = c("up", "down", "un")
)

ggplot(nData, aes(x = S2D16_pro, y = ADSC_pro)) +
  geom_tile(aes(fill = type), color = "black", size = .5, show.legend = F) +
  geom_text(aes(label = N)) +
  scale_fill_manual(values = S2col) +
  coord_fixed() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("pro_N.pdf", width = 3, height = 3)

S2pro <- allPeakMeta[
  type == "pro" & (H3K4me3_group != "-" | H3K27me3_group != "-") &
    S2D16_pro == "Act" & ADSC_pro %ni% c("Act"), id]

ADSCpro <- allPeakMeta[
  type == "pro" & (H3K4me3_group != "-" | H3K27me3_group != "-") &
    ADSC_pro == "Act" & S2D16_pro %ni% c("Act"), id]

geneOver <- allPeakMeta[id %in% c(S2pro, ADSCpro)]
geneOver$type <- "ADSC"
geneOver[id %in% S2pro, type := "S2"][]
table(geneOver$type)

geneList <- list(
  S2 = geneOver[type == "S2", gene] %>% str_c(collapse = "|") %>% str_split("\\|") %>% pluck(1) %>% unique(),
  ADSC = geneOver[type == "ADSC", gene] %>% str_c(collapse = "|") %>% str_split("\\|") %>% pluck(1) %>% unique()
)

saveRDS(geneList, "diffPro_geneList.rds")

# gene list signal --------------------------------------------------------

fsdir <- "dediff"
fs <- readRDS("fs.rds")

proData <- fread("geneBody.bed")
proData[, V1 := str_remove(V1, "chr")]
fsPro <- proData[V4 %in% fs]

fwrite(fsPro, glue("{fsdir}/dediff.bed"), sep = "\t", col.names = F)

gbwDir <- "zbw"

deepDir <- "deeptools_env/bin"
fileName <- c()
fileName <- list.files(gbwDir, ".bw", full.names = T) %>% c(fileName)

usedHis <- c("H3K4me3", "H3K27me3")
fileList <- map(usedHis, ~ str_subset(fileName, .x)) %>% set_names(usedHis)

proCMD <- c()
for(i in names(fileList)) {
  message(i)
  proCMD <- glue(
    "{deepDir}/computeMatrix reference-point -S {bws} -R dediff.bed \\
    -a 2000 -b 2000 -bs 20 --referencePoint TSS -p 12 \\
    -o {i}.mat --outFileNameMatrix {i}.tsv --outFileSortedRegions {i}.bed \\\n
    
    {deepDir}/plotHeatmap -m {i}.mat --colorMap RdBu_r \\
    -o {i}.heatmap.pdf --outFileSortedRegions {i}.sort.bed",
    bws = fileList[[i]] %>% str_c(collapse = " ")) %>% c(proCMD, .)
}
cat(proCMD[1])

write.table(c("#!/bin/bash", proCMD), glue("{fsdir}/pro.sh"), sep = "\n", row.names = F, col.names = F, quote = F)

for(i in c("H3K4me3", "H3K27me3")) {
  message(i)
  markData <- fread(str_c(fsdir, "/", i, ".tsv"), skip = 3, header = F) %>% as.matrix()
  markMeta <- fread(str_c(fsdir, "/", i, ".bed"))
  rownames(markData) <- markMeta$name
  markData[is.na(markData)] <- 0
  
  proMean <- markData %>% apply(1, function(x) tapply(x, rep(fileList[[i]] %>% str_remove("_zscore.*"), each = 200), mean)) %>% t
  rn <- colnames(proMean) %>% str_remove("_H.*")
  proMean <- apply(proMean, 1, function(x) tapply(x, rn, mean)) %>% t
  saveRDS(proMean, str_c(fsdir, "/", i, "_proMean.rds"))
}

# after 'enhancer' part of scATAC+RNA_hADSC_GRN.R

gs <- "dediff"
gs <- "fib"

gl <- readRDS("S2key.rds")
gl <- readRDS("ADSCkey.rds")

e2gList <- readRDS("e2gList.rds")

gl <- intersect(gl, names(e2gList))

markData <- list()

markData$H3K4me3 <- readRDS(str_c(gs, "/H3K4me3_proMean.rds"))
markData$H3K27me3 <- readRDS(str_c(gs, "/H3K27me3_proMean.rds"))
markData$H3K4me3[markData$H3K4me3 < 0] <- 0
markData$H3K27me3[markData$H3K27me3 < 0] <- 0

e2g <- e2gList[gl]
zGene <- gl[map_lgl(e2g, ~ identical(.x, "none"))]
e2g <- e2g[setdiff(gl, zGene)]

e2g <- imap(e2g, ~ .x[, gene := .y])
e2g %<>% do.call(rbind, .)

e2g[, w := cor/mean(cor), by = gene]

i = "H3K27me3"
cn <- colnames(valList$H3K4me3) %>% str_remove_all("_H.*")

enhData <- list()

for(i in c("H3K4me1", "H3K27ac")) {
  message(i)
  peakList <- allPeakMeta[match(e2g$peak, id), ..i][[1]] %>% str_split(";")
  enhData[[i]] <- map(peakList, ~ {
    if(identical(.x, "-")) {rep(0, length(unique(cn)))} else
      if(length(.x) == 1) {
        l <- peakL(.x)
        w <- l / 1000
        val <- tapply(valList[[i]][.x, ], cn, mean)[unique(cn)]
        val[val < 0] <- 0
        val * w
      } else {
        l <- peakL(.x)
        w <- l / 1000
        valMtx <- apply(valList[[i]][.x, ], 1, function(x) tapply(x, cn, mean)) %>% t
        valMtx[valMtx < 0] <- 0
        valMtx <- valMtx[, unique(cn)]
        sweep(valMtx, 1, w, "*") %>% colMeans
      }
  }) %>% do.call(rbind, .)
}

enhData$H3K4me1 <- enhData$H3K4me1 * e2g$w
enhData$H3K27ac <- enhData$H3K27ac * e2g$w

markData$H3K4me1 <- apply(enhData$H3K4me1, 2, function(x) tapply(x, e2g$gene, sum))
markData$H3K27ac <- apply(enhData$H3K27ac, 2, function(x) tapply(x, e2g$gene, sum))

for(i in c("H3K4me1", "H3K27ac")) {
  noMat <- matrix(0, nrow = length(zGene), ncol = ncol(markData[[i]])) %>% set_rownames(zGene)
  markData[[i]] <- rbind(noMat, markData[[i]])
}

markData <- map(markData, ~ .x[gl, ])

saveRDS(markData, glue("{gs}/markData.rds"))
markData <- readRDS(glue("{gs}/markData.rds"))

th <- readRDS("th.rds")

# enh K27me3
markData <- list()
enhData <- list()

for(i in c("H3K27me3")) {
  message(i)
  peakList <- allPeakMeta[match(e2g$peak, id), ..i][[1]] %>% str_split(";")
  enhData[[i]] <- map(peakList, ~ {
    if(identical(.x, "-")) {rep(0, length(unique(cn)))} else
      if(length(.x) == 1) {
        l <- peakL(.x)
        w <- l / 1000
        val <- tapply(valList[[i]][.x, ], cn, mean)[unique(cn)]
        val[val < 0] <- 0
        val * w
      } else {
        l <- peakL(.x)
        w <- l / 1000
        valMtx <- apply(valList[[i]][.x, ], 1, function(x) tapply(x, cn, mean)) %>% t
        valMtx[valMtx < 0] <- 0
        valMtx <- valMtx[, unique(cn)]
        sweep(valMtx, 1, w, "*") %>% colMeans
      }
  }) %>% do.call(rbind, .)
}

enhData$H3K27me3 <- enhData$H3K27me3 * e2g$w
markData$H3K27me3 <- apply(enhData$H3K27me3, 2, function(x) tapply(x, e2g$gene, sum))

for(i in "H3K27me3") {
  noMat <- matrix(0, nrow = length(zGene), ncol = ncol(markData[[i]])) %>% set_rownames(zGene)
  markData[[i]] <- rbind(noMat, markData[[i]])
}
markData <- map(markData, ~ .x[gl, ])

saveRDS(markData, glue("{gs}/markData_enhK27me3.rds"))

# load
markData <- readRDS(glue("{gs}/markData.rds"))
markData2 <- readRDS(glue("{gs}/markData_enhK27me3.rds"))
names(markData2) <- "H3K27me3_enh"

markData <- c(markData, markData2)

markCol <- structure(
  c(getMatColor("red", 6)[6],
    getMatColor("orange", 4)[4],
    getMatColor("teal", 1),
    getMatColor("indigo", 1),
    getMatColor("blue", 1)),
  names = c("H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3", "H3K27me3_enh")
)
markCol %>% show_col()

lfcData <- list()
lfcAvg <- list()
lfcSe <- list()

for(i in names(markData)) {
  lfc <- log2(markData[[i]] + 1)
  colnames(lfc) %<>% str_remove(".*/")
  lfc <- lfc - lfc[, "ADSC"]
  lfc <- lfc[, sort(colnames(lfc))]
  lfcData[[i]] <- lfc
  lfcAvg[[i]] <- colMeans(lfc)
  lfcSe[[i]] <- apply(lfc, 2, function(x) sd(x)/sqrt(length(x)))
}

plotData <- do.call(cbind, lfcAvg) %>% reshape2::melt() %>% set_colnames(c("samp", "mark", "value")) %>% as.data.table
plotData$se <- do.call(cbind, lfcSe) %>% reshape2::melt() %>% set_colnames(c("samp", "mark", "se")) %>% {.[, "se"]}

samp <- c("ADSC", "S1D4", "S1D8", "S2D4", "S2D8", "S2D12", "S2D16", "IPSno3", "IPSno7", "H1", "H9")
plotData$samp %<>% factor(samp)
plotData <- plotData[order(samp)]

plotData[, x := as.numeric(samp)][]
plotData[mark == "H3K4me3", x := x - 0.1][]
plotData[mark == "H3K27ac", x := x - 0.05][]
plotData[mark == "H3K4me1", x := x + 0][]
plotData[mark == "H3K27me3", x := x + 0.05][]
plotData[mark == "H3K27me3_enh", x := x + 0.1][]

plotData$mark %<>% factor(names(markCol))
plotData <- plotData[order(mark)]

lineData <- plotData[samp %ni% c("IPSno3", "IPSno7", "H1", "H9")]

ggplot(plotData, aes(x = x, y = value, group = mark)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_path(
  #   aes(color = mark), data = lineData,
  #   alpha = .5, show.legend = F) +
  geom_point(aes(color = mark), shape = 15, size = 1) +
  geom_linerange(aes(color = mark, ymin = value - se, ymax = value + se), show.legend = F) +
  scale_color_manual(values = markCol) +
  scale_x_continuous(
    breaks = seq_along(samp),
    labels = samp) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "", y = "Log2FC to ADSC", title = gs) +
  # labs(x = "", y = "Relative signal to ADSC", title = gs) +
  theme(
    aspect.ratio = 1/2,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    legend.key = element_blank()
  )

ggsave(str_c(gs, "_mark_lfc.pdf"), width = 7, height = 3)

# gene set typing ---------------------------------------------------------

# after 'enhancer' part of scATAC+RNA_hADSC_GRN.R

e2gList <- readRDS("e2gList.rds")

gs <- "dediff"
gl <- readRDS("S2key.rds")
gs <- "fib"
gl <- readRDS("ADSCkey.rds")

markData <- list()

markData$H3K4me3 <- readRDS(str_c(gs, "/H3K4me3_proMean.rds"))
markData$H3K27me3 <- readRDS(str_c(gs, "/H3K27me3_proMean.rds"))
markData$H3K4me3[markData$H3K4me3 < 0] <- 0
markData$H3K27me3[markData$H3K27me3 < 0] <- 0

e2g <- e2gList[gl]
zGene <- gl[map_lgl(e2g, ~ identical(.x, "none"))]
e2g <- e2g[setdiff(gl, zGene)]

e2g <- imap(e2g, ~ .x[, gene := .y])
e2g %<>% do.call(rbind, .)

e2g[, w := cor/mean(cor), by = gene]

i = "H3K27ac"
cn <- colnames(valList$H3K4me3) %>% str_remove_all("_H.*")

enhData <- list()

for(i in c("H3K4me1", "H3K27ac")) {
  message(i)
  peakList <- allPeakMeta[match(e2g$peak, id), ..i][[1]] %>% str_split(";")
  enhData[[i]] <- map(peakList, ~ {
    if(identical(.x, "-")) {rep(0, length(unique(cn)))} else
      if(length(.x) == 1) {
        l <- peakL(.x)
        w <- l / 1000
        val <- tapply(valList[[i]][.x, ], cn, mean)[unique(cn)]
        val[val < 0] <- 0
        val * w
      } else {
        l <- peakL(.x)
        w <- l / 1000
        valMtx <- apply(valList[[i]][.x, ], 1, function(x) tapply(x, cn, mean)) %>% t
        valMtx[valMtx < 0] <- 0
        valMtx <- valMtx[, unique(cn)]
        sweep(valMtx, 1, w, "*") %>% colMeans
      }
  }) %>% do.call(rbind, .)
}

enhData$H3K4me1 <- enhData$H3K4me1 * e2g$w
enhData$H3K27ac <- enhData$H3K27ac * e2g$w

markData$H3K4me1 <- apply(enhData$H3K4me1, 2, function(x) tapply(x, e2g$gene, sum))
markData$H3K27ac <- apply(enhData$H3K27ac, 2, function(x) tapply(x, e2g$gene, sum))

for(i in c("H3K4me1", "H3K27ac")) {
  noMat <- matrix(0, nrow = length(zGene), ncol = ncol(markData[[i]])) %>% set_rownames(zGene)
  markData[[i]] <- rbind(noMat, markData[[i]])
}

markData <- map(markData, ~ .x[gl, ])

saveRDS(markData, glue("{gs}/markData.rds"))
markData <- readRDS(glue("{gs}/markData.rds"))

th <- readRDS("th.rds")

# pro
glPro <- data.table(
  gene = gl,
  K27me3_S2 = log2(markData$H3K27me3[, "S2D16"] + 1) - log2(markData$H3K27me3[, "ADSC"] + 1),
  K27me3_S1 = log2(markData$H3K27me3[, "S1D8"] + 1) - log2(markData$H3K27me3[, "ADSC"] + 1),
  K4me3_S2 = log2(markData$H3K4me3[, "S2D16"] + 1) - log2(markData$H3K4me3[, "ADSC"] + 1),
  K4me3_S1 = log2(markData$H3K4me3[, "S1D8"] + 1) - log2(markData$H3K4me3[, "ADSC"] + 1)
)

for(k in c("S2D16", "S1D8", "ADSC")) {
  k4 <- markData$H3K4me3[glPro$gene, k] > th$H3K4me3
  k27 <- markData$H3K27me3[glPro$gene, k] > th$H3K27me3
  set(glPro, j = k, value = "NoSig")
  set(glPro, i = which(k4 & k27), j = k, value = "Biv")
  set(glPro, i = which(k4 & !k27), j = k, value = "Act")
  set(glPro, i = which(!k4 & k27), j = k, value = "Ina")
}

table(glPro[, c("S1D8", "ADSC")])
table(glPro[, c("S2D16", "ADSC")])
ggplot(glPro, aes(x = K4me3_S1, y = K27me3_S1)) +
  geom_point()
ggplot(glPro, aes(x = K4me3_S2, y = K27me3_S2)) +
  geom_point()
ggplot(glPro, aes(x = K4me3_S1, y = K4me3_S2)) +
  geom_point()

glPro[gene %in% c("SALL4", "LIN28A", "HOXB9", "MSX1", "MSX2", "MYCN")]


# type S2 (pro)
glPro$S2type <- "other"
glPro[ADSC %in% c("NoSig", "Ina") & S2D16 %in% c("NoSig", "Ina"), S2type := "noK4"]
glPro[ADSC %in% c("Act") & S2D16 %in% c("Act"), S2type := "noK27"]

table(glPro$S2type)

ggplot(glPro, aes(x = K4me3_S2, y = K27me3_S2)) +
  geom_point(aes(color = S2type))

K4avg <- glPro[S2type == "noK4", K4me3_S2] %>% mean
K4sd <- glPro[S2type == "noK4", K4me3_S2] %>% sd
K27avg <- glPro[S2type == "noK27", K27me3_S2] %>% mean
K27sd <- glPro[S2type == "noK27", K27me3_S2] %>% sd

glPro$group <- "other"
glPro[K27me3_S2 < K27avg - K27sd & K4me3_S2 > K4avg + K4sd, group := "core"]
glPro[, score := K4me3_S2/sd(K4me3_S2) - K27me3_S2/sd(K27me3_S2)]
glPro[order(-score)]

glPro[group == "core"][order(-K4me3_S2)]

write.csv(glPro, glue("{gs}_type.csv"), row.names = F)
saveRDS(glPro, "dediff_glPro.rds")
glPro <- readRDS("dediff_glPro.rds")

xRange <- range(glPro$K4me3_S2)
yRange <- range(glPro$K27me3_S2)
xmin <- xRange[1] - diff(xRange) * .05
xmax <- xRange[2] + diff(xRange) * .05
ymin <- yRange[1] - diff(yRange) * .05
ymax <- yRange[2] + diff(yRange) * .05

groupCol <- c(
  core = getMatColor("red", 5)[5],
  other = "gray60"
)

pctData <- glPro[, .N, by = group]
pct <- pctData[, round(N/sum(N) * 100, 2)][pctData$group == "core"] %>% str_c("%")

ggplot(glPro[order(score)], aes(x = K4me3_S2, y = K27me3_S2)) +
  geom_rect(xmin = K4avg - K4sd, xmax = K4avg + K4sd, ymin = ymin, ymax = ymax, fill = "gray90", color = NA) +
  geom_rect(ymin = K27avg - K27sd, ymax = K27avg + K27sd, xmin = xmin, xmax = xmax, fill = "gray90", color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point_rast(aes(color = group), shape = 16, size = .5) +
  geom_text_repel(
    data = glPro[gene %in% c("SALL4", "LIN28A", "HOXB9", "MYCN")],
    aes(label = gene), size = 2, fontface = "italic",
    color = "black", point.padding = .1,
    segment.size = 0.2, min.segment.length = 0.1,
    nudge_x = .5, nudge_y = -0.5) +
  scale_color_manual(values = groupCol) +
  scale_x_continuous(limits = c(xmin, xmax), expand = expansion(c(0, 0))) +
  scale_y_continuous(limits = c(ymin, ymax), expand = expansion(c(0, 0))) +
  labs(x = "H3K4me3 Log2FC", y = "H3K27me3 Log2FC", color = "", title = "S2D16.vs.ADSC", tag = pct) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    plot.tag = element_text(size = 8, color = groupCol[1]),
    plot.tag.position = c(0.7, 0.35),
    legend.key = element_blank()
  )
ggsave(glue("{gs}_proDiff.pdf"), width = 6, height = 4, scale = .7)

# type S1 (pro)
glPro$S1type <- "other"
glPro[ADSC %in% c("NoSig", "Ina") & S1D8 %in% c("NoSig", "Ina"), S1type := "noK4"]
glPro[ADSC %in% c("Act") & S1D8 %in% c("Act"), S1type := "noK27"]

table(glPro$S1type)

ggplot(glPro, aes(x = K4me3_S1, y = K27me3_S1)) +
  geom_point(aes(color = S1type))

K4avg <- glPro[S1type == "noK4", K4me3_S1] %>% mean
K4sd <- glPro[S1type == "noK4", K4me3_S1] %>% sd
K27avg <- glPro[S1type == "noK27", K27me3_S1] %>% mean
K27sd <- glPro[S1type == "noK27", K27me3_S1] %>% sd

glPro$S1group <- "other"
glPro[K27me3_S1 < K27avg - K27sd & K4me3_S1 > K4avg + K4sd, S1group := "core"]
glPro[, S1score := K4me3_S1/sd(K4me3_S1) - K27me3_S1/sd(K27me3_S1)]
glPro[order(-S1score)]

glPro[S1group == "core"][order(-K4me3_S1)]

xRange <- range(glPro$K4me3_S1)
yRange <- range(glPro$K27me3_S1)
xmin <- xRange[1] - diff(xRange) * .05
xmax <- xRange[2] + diff(xRange) * .05
ymin <- yRange[1] - diff(yRange) * .05
ymax <- yRange[2] + diff(yRange) * .05

groupCol <- c(
  core = getMatColor("red", 5)[2],
  other = "gray60"
)
groupCol %>% show_col()

pctData <- glPro[, .N, by = S1group]
pct <- pctData[, round(N/sum(N) * 100, 2)][pctData$S1group == "core"] %>% str_c("%")

ggplot(glPro[order(S1score)], aes(x = K4me3_S1, y = K27me3_S1)) +
  geom_rect(xmin = K4avg - K4sd, xmax = K4avg + K4sd, ymin = ymin, ymax = ymax, fill = "gray90", color = NA) +
  geom_rect(ymin = K27avg - K27sd, ymax = K27avg + K27sd, xmin = xmin, xmax = xmax, fill = "gray90", color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point_rast(aes(color = S1group), shape = 16, size = .5) +
  geom_text_repel(
    data = glPro[gene %in% c("SALL4", "LIN28A", "HOXB9", "MYCN")],
    aes(label = gene), size = 2, fontface = "italic",
    color = "black", point.padding = .1,
    segment.size = 0.2, min.segment.length = 0.1,
    nudge_x = .5, nudge_y = -0.3) +
  scale_color_manual(values = groupCol) +
  scale_x_continuous(limits = c(xmin, xmax), expand = expansion(c(0, 0))) +
  scale_y_continuous(limits = c(ymin, ymax), expand = expansion(c(0, 0))) +
  labs(x = "H3K4me3 Log2FC", y = "H3K27me3 Log2FC", color = "", title = "S1D8.vs.ADSC", tag = pct) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    plot.tag = element_text(size = 8, color = groupCol[1]),
    plot.tag.position = c(0.7, 0.4),
    legend.key = element_blank()
  )
ggsave(glue("{gs}_proDiff_S1.pdf"), width = 6, height = 4, scale = .7)


# enh
glEnh <- data.table(
  gene = gl,
  K27ac_S2 = log2(markData$H3K27ac[, "S2D16"] + 1) - log2(markData$H3K27ac[, "ADSC"] + 1),
  K27ac_S1 = log2(markData$H3K27ac[, "S1D8"] + 1) - log2(markData$H3K27ac[, "ADSC"] + 1),
  K4me1_S2 = log2(markData$H3K4me1[, "S2D16"] + 1) - log2(markData$H3K4me1[, "ADSC"] + 1),
  K4me1_S1 = log2(markData$H3K4me1[, "S1D8"] + 1) - log2(markData$H3K4me1[, "ADSC"] + 1)
)

for(k in c("S2D16", "S1D8", "ADSC")) {
  k4 <- markData$H3K4me1[glEnh$gene, k] > th$H3K4me1
  k27 <- markData$H3K27ac[glEnh$gene, k] > th$H3K27ac
  set(glEnh, j = k, value = "NoSig")
  set(glEnh, i = which(k4 & k27), j = k, value = "Act")
  set(glEnh, i = which(k4 & !k27), j = k, value = "Poi")
  set(glEnh, i = which(!k4 & k27), j = k, value = "Nor")
}

table(glEnh[, c("S1D8", "ADSC")])
table(glEnh[, c("S2D16", "ADSC")])
ggplot(glEnh, aes(x = K4me1_S1, y = K27ac_S1)) +
  geom_point()
ggplot(glEnh, aes(x = K4me1_S2, y = K27ac_S2)) +
  geom_point()
ggplot(glEnh, aes(x = K4me1_S1, y = K4me1_S2)) +
  geom_point()

glEnh[gene %in% c("SALL4", "LIN28A", "HOXB9", "MSX1", "MSX2", "MYCN")]

# type S2 (enh)
glEnh$S2type <- "other"
glEnh[ADSC %in% c("NoSig") & S2D16 %in% c("NoSig"), S2type := "noSig"]

table(glEnh$S2type)

ggplot(glEnh, aes(x = K4me1_S2, y = K27ac_S2)) +
  geom_point(aes(color = S2type))

K4avg <- glEnh[S2type == "noSig", K4me1_S2] %>% mean
K4sd <- glEnh[S2type == "noSig", K4me1_S2] %>% sd
K27avg <- glEnh[S2type == "noSig", K27ac_S2] %>% mean
K27sd <- glEnh[S2type == "noSig", K27ac_S2] %>% sd

glEnh$group <- "other"
glEnh[K27ac_S2 > K27avg + K27sd & K4me1_S2 > K4avg + K4sd, group := "core"]
glEnh[, score := K4me1_S2/sd(K4me1_S2) + K27ac_S2/sd(K27ac_S2)]
glEnh[order(-score)]

glEnh[group == "core"][order(-K4me1_S2)]

write.csv(glEnh, glue("{gs}_Enh_type.csv"), row.names = F)
saveRDS(glEnh, "dediff_glEnh.rds")
glEnh <- readRDS("dediff_glEnh.rds")

saveRDS(glEnh, "mj_glEnh.rds")

xRange <- range(glEnh$K4me1_S2)
yRange <- range(glEnh$K27ac_S2)
xmin <- xRange[1] - diff(xRange) * .05
xmax <- xRange[2] + diff(xRange) * .05
ymin <- yRange[1] - diff(yRange) * .05
ymax <- yRange[2] + diff(yRange) * .05

groupCol <- c(
  core = getMatColor("blue", 5)[5],
  other = "gray60"
)

pctData <- glEnh[, .N, by = group]
pct <- pctData[, round(N/sum(N) * 100, 2)][pctData$group == "core"] %>% str_c("%")

ggplot(glEnh[order(score)], aes(x = K4me1_S2, y = K27ac_S2)) +
  geom_rect(xmin = K4avg - K4sd, xmax = K4avg + K4sd, ymin = ymin, ymax = ymax, fill = "gray90", color = NA) +
  geom_rect(ymin = K27avg - K27sd, ymax = K27avg + K27sd, xmin = xmin, xmax = xmax, fill = "gray90", color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point_rast(aes(color = group), shape = 16, size = .5) +
  geom_text_repel(
    data = glEnh[gene %in% c("SALL4", "LIN28A", "MSX1", "MSX2", "MYCN")],
    aes(label = gene), size = 2, fontface = "italic",
    color = "black", point.padding = .1,
    segment.size = 0.2, min.segment.length = 0.1,
    nudge_x = .5, nudge_y = -0.5) +
  scale_color_manual(values = groupCol) +
  scale_x_continuous(limits = c(xmin, xmax), expand = expansion(c(0, 0))) +
  scale_y_continuous(limits = c(ymin, ymax), expand = expansion(c(0, 0))) +
  labs(x = "H3K4me1 Log2FC", y = "H3K27ac Log2FC", color = "", title = "S2D16.vs.ADSC", tag = pct) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    plot.tag = element_text(size = 8, color = groupCol[1]),
    plot.tag.position = c(0.7, 0.65),
    legend.key = element_blank()
  )
ggsave(glue("{gs}_enhDiff.pdf"), width = 6, height = 4, scale = .7)

# type S1 (enh)
glEnh$S1type <- "other"
glEnh[ADSC %in% c("NoSig") & S1D8 %in% c("NoSig"), S1type := "noSig"]

table(glEnh$S1type)

ggplot(glEnh, aes(x = K4me1_S1, y = K27ac_S1)) +
  geom_point(aes(color = S1type))

K4avg <- glEnh[S1type == "noSig", K4me1_S1] %>% mean
K4sd <- glEnh[S1type == "noSig", K4me1_S1] %>% sd
K27avg <- glEnh[S1type == "noSig", K27ac_S1] %>% mean
K27sd <- glEnh[S1type == "noSig", K27ac_S1] %>% sd

glEnh$S1group <- "other"
glEnh[K27ac_S1 > K27avg + K27sd & K4me1_S1 > K4avg + K4sd, S1group := "core"]
glEnh[, S1score := K4me1_S1/sd(K4me1_S1) + K27ac_S1/sd(K27ac_S1)]
glEnh[order(-S1score)]

glEnh[S1group == "core"][order(-K4me1_S1)]

xRange <- range(glEnh$K4me1_S1)
yRange <- range(glEnh$K27ac_S1)
xmin <- xRange[1] - diff(xRange) * .05
xmax <- xRange[2] + diff(xRange) * .05
ymin <- yRange[1] - diff(yRange) * .05
ymax <- yRange[2] + diff(yRange) * .05

groupCol <- c(
  core = getMatColor("blue", 5)[2],
  other = "gray60"
)
groupCol %>% show_col()

pctData <- glEnh[, .N, by = S1group]
pct <- pctData[, round(N/sum(N) * 100, 2)][pctData$S1group == "core"] %>% str_c("%")

ggplot(glEnh[order(score)], aes(x = K4me1_S1, y = K27ac_S1)) +
  geom_rect(xmin = K4avg - K4sd, xmax = K4avg + K4sd, ymin = ymin, ymax = ymax, fill = "gray90", color = NA) +
  geom_rect(ymin = K27avg - K27sd, ymax = K27avg + K27sd, xmin = xmin, xmax = xmax, fill = "gray90", color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point_rast(aes(color = S1group), shape = 16, size = .5) +
  geom_text_repel(
    data = glEnh[gene %in% c("SALL4", "LIN28A", "MSX1", "MSX2", "MYCN")],
    aes(label = gene), size = 2, fontface = "italic",
    color = "black", point.padding = .1,
    segment.size = 0.2, min.segment.length = 0.1,
    nudge_x = .75, nudge_y = -0.5) +
  scale_color_manual(values = groupCol) +
  scale_x_continuous(limits = c(xmin, xmax), expand = expansion(c(0, 0))) +
  scale_y_continuous(limits = c(ymin, ymax), expand = expansion(c(0, 0))) +
  labs(x = "H3K4me1 Log2FC", y = "H3K27ac Log2FC", color = "", title = "S1D8.vs.ADSC", tag = pct) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    plot.tag = element_text(size = 8, color = groupCol[1]),
    plot.tag.position = c(0.7, 0.75),
    legend.key = element_blank()
  )
ggsave(glue("{gs}_enhDiff_S1.pdf"), width = 6, height = 4, scale = .7)

# rect venn plot ----------------------------------------------------------

table(glPro[, c("group")])
table(glPro[, c("S1group", "group")])
table(glEnh[, c("S1group", "group")])

identical(glEnh$gene, glPro$gene)
glEnh$proS2 <- glPro$group
table(glEnh[, c("proS2", "S1group")])
table(glEnh[, c("proS2", "group")])

usedGene <- intersect(glPro[S1group == "other", gene], glEnh[S1group == "core", gene])
glEnh[gene %in% usedGene & proS2 == "core"][order(-K27ac_S1)]

glEnh[order(-K27ac_S1)]

plotRectVenn(
  pct = 33/77,
  txt = c(44, 3, 33),
  label = c("S2 pro", "S1 pro"),
  typeCol = c(
    getMatColor("red", 5)[c(5, 2)],
    getMatColor("red", 5)[3]
  )
)
plotRectVenn(
  pct = (32+23+88)/(32+23+88+18),
  txt = c(18, 23, 32+23+88),
  label = c("S2 pro", "S1 pro"),
  typeCol = c(
    getMatColor("indigo", 5)[c(3)],
    getMatColor("cyan", 5)[3],
    getMatColor("light-blue", 5)[3]
  )
)
ggsave("fib_S2pro_S1pro.venn.pdf", width = 5, height = 5)
plotRectVenn(
  pct = 81/(81+13),
  txt = c(13, 11, 81),
  label = c("S2 enh", "S1 enh"),
  typeCol = c(
    getMatColor("blue", 5)[c(5, 2, 3)]
  )
)
plotRectVenn(
  pct = 127/(127+25),
  txt = c(25, 2, 127),
  label = c("S2 enh", "S1 enh"),
  typeCol = c(
    getMatColor("light-green", 5)[c(5, 2, 3)]
  )
)
ggsave("fib_S2enh_S1enh.venn.pdf", width = 5, height = 5)

plotRectVenn(
  pct = 62/(62+15),
  txt = c(15, 30, 62),
  label = c("S2 pro", "S1 enh"),
  typeCol = c(
    getMatColor("red", 5)[c(5)],
    getMatColor("blue", 5)[c(2)],
    getMatColor("purple", 5)[c(2)]
  )
)
plotRectVenn(
  pct = 57/(57+20),
  txt = c(20, 37, 57),
  label = c("S2 pro", "S2 enh"),
  typeCol = c(
    getMatColor("red", 5)[c(5)],
    getMatColor("blue", 5)[c(5)],
    getMatColor("purple", 5)[c(4)]
  )
)

ggsave("S2pro_S1pro.venn.pdf", width = 5, height = 5)
ggsave("S2enh_S1enh.venn.pdf", width = 5, height = 5)
ggsave("S2pro_S1enh.venn.pdf", width = 5, height = 5)
ggsave("S2pro_S2enh.venn.pdf", width = 5, height = 5)

# track plot --------------------------------------------------------------

library(rtracklayer)
library(patchwork)
library(gtable)

transGene <- readRDS("hg19canonicalTrans.rds")
transMeta <- readRDS("hg19trans.rds") %>% as.data.table()
exonMeta <- readRDS("hg19exon.rds") %>% as.data.table()

bwFile <- list.files("mzbw", recursive = T, full.names = T)

usedHis <- c("H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3")
sel <- map(usedHis, ~ str_subset(bwFile, .x)) %>% set_names(usedHis)

e2gList <- readRDS("e2gList.rds")
cA <- readRDS("ao_cA.rds")
allPeakMeta <- readRDS("allPeakMeta.rds")

genePro <- fread("promoter2k_sort.bed")
genePro %<>% set_colnames(c("chr", "start", "end", "gene"))

sn <- c(
  "ADSC",
  "S1D4", "S1D8",
  "S2D4", "S2D8", "S2D12", "S2D16",
  "IPSno3", "IPSno7", "H1", "H9")
sampOrder <- c(3, 6:7, 10:11, 8:9, 5, 4, 1:2)

markCol <- structure(
  c(getMatColor("red", 6)[6],
    getMatColor("orange", 4)[4],
    getMatColor("teal", 1),
    getMatColor("indigo", 1)),
  names = c("H3K4me3", "H3K27ac", "H3K4me1", "H3K27me3")
)
markCol %>% show_col()

sampColor <- structure(c(
  getMatColor("blue", 1),
  getMatColor("green", 2),
  getMatColor("yellow", 1),
  getMatColor("amber", 1),
  getMatColor("orange", 1),
  getMatColor("deep-orange", 1),
  getMatColor("red", 2),
  getMatColor("purple", 2)), names = c(
    "ADSC", "S1D4", "S1D8",
    "S2D4", "S2D8", "S2D12", "S2D16",
    "IPSno3", "IPSno7", "H1", "H9"))
sampColor %>% show_col()

hasEnh <- map_lgl(e2gList, ~ identical(.x, "none"))
noGene <- names(hasEnh)[hasEnh]

usedGene <- c("")

figList <- list()

for (gene in usedGene) {
  message(gene)
  
  gAnno <- getGeneAnno(gene)
  canPlotLink <- gene %ni% noGene
  
  if(!canPlotLink) {
    usedPeak <- character(0)
  } else {usedPeak <- e2gList[[gene]]$peak}
  selReg <- allPeakMeta[id %in% usedPeak, 1:4]
  selReg$color <- "gray90"
    
  isG <- genePro$gene == gene
  selReg <- cbind(genePro[isG], color = getMatColor("yellow", 2)[1]) %>% unname %>% rbind(selReg)
  
  gLink <- cA[seqnames == str_c("chr", gAnno$geneChr)]
  gLink[, keep := "no"]
  
  for(i in 1:nrow(selReg)) {
    gLink[(start >= selReg[i, start] & start <= selReg[i, end]) | (end >= selReg[i, start] & end <= selReg[i, end]), keep := "yes"]
  }
  gLink <- gLink[keep == "yes"]
  
  pList <- list()
  for(i in names(sel)) {
    # message(i)
    pList[[i]] <- plotTrack(
      chr = gAnno$geneChr,
      startPos = gAnno$startPos,
      endPos = gAnno$endPos,
      extendLength = 2e3,
      nbin = 2000,
      bwFile = sel[[i]][sampOrder],
      sampleName = sn,
      sampleColor = markCol[i] %>% unname %>% rep(length(sn)),
      stripColor = sampColor,
      addRegion = T, selRegion = selReg,
      showName = F, showTh = F,
      ratio = 1/10)
  }
  
  pList$anno <- plotAnno(
    chr = gAnno$geneChr,
    startPos = min(c(gAnno$startPos, selReg$start)),
    endPos = max(c(gAnno$endPos, selReg$end)),
    extendLength = 2e3,
    gene = gene,
    ratio = 0.15, dpi = 600
  )
  
  if(canPlotLink) {
    pList$link <- plotLink(
      chr = gAnno$geneChr, addChr = T,
      startPos = min(c(gAnno$startPos, selReg$start)),
      endPos = max(c(gAnno$endPos, selReg$end)),
      extendLength = 2e3,
      nbin = 2000,
      linkData = gLink,
      ratio = 0.15, dpi = 600
    )
  }
  
  for(i in names(sel)) {
    pList[[i]]$widths[pList[[i]]$layout$l[str_detect(pList[[i]]$layout$name, "axis-l")] %>% unique] <- unit(0, "null")
    pList[[i]]$heights[pList[[i]]$layout$t[str_detect(pList[[i]]$layout$name, "axis-b")]] <- unit(0, "null")
    pList[[i]]$widths[pList[[i]]$layout$l[str_detect(pList[[i]]$layout$name, "strip-r")] %>% unique] <- unit(0.04, "null")
    pList[[i]]$widths[length(pList[[i]]$widths)] <- unit(0.02, "null")
  }
  
  pList$anno <- ggplotGrob(pList$anno)
  pList$anno$widths[pList$anno$layout$l[str_detect(pList$anno$layout$name, "axis-l")] %>% unique] <- unit(0, "null")
  pList$anno$widths[pList$anno$layout$l[str_detect(pList$anno$layout$name, "ylab-l")] %>% unique] <- unit(0, "null")
  pList$anno$heights[pList$anno$layout$t[str_detect(pList$anno$layout$name, "axis-b")] %>% unique] <- unit(0, "null")
  pList$anno$heights[pList$anno$layout$t[str_detect(pList$anno$layout$name, "xlab-b")] %>% unique] <- unit(0, "null")
  pList$anno$widths[length(pList$anno$widths)] <- unit(0.06, "null")
  
  if(canPlotLink) {
    pList$link <- ggplotGrob(pList$link)
    pList$link$widths[pList$link$layout$l[str_detect(pList$link$layout$name, "axis-l")] %>% unique] <- unit(0, "null")
    pList$link$heights[pList$link$layout$t[str_detect(pList$link$layout$name, "axis-b")] %>% unique] <- unit(0, "null")
    pList$link$widths[length(pList$link$widths)] <- unit(0.06, "null")
  }
  
  fig <- gtable(
    widths = unit(rep(1.06, length(sel)), "null"),
    heights = unit(c(0.1 * length(sn), 0.15, 0.15), "null"))
  
  for(i in 1:length(sel)) {
    fig <- gtable_add_grob(
      fig,
      grobs = pList[[i]],
      t = 1,
      l = i)
    
    if(canPlotLink) {
      fig <- gtable_add_grob(
        fig,
        grobs = pList$link,
        t = 2,
        l = i)
    }
    
    fig <- gtable_add_grob(
      fig,
      grobs = pList$anno,
      t = 3,
      l = i)
  }
  
  figList[[gene]] <- fig
}

grid.newpage()
grid.draw(figList$SALL4)
grid.draw(figList$LEF1)

for(gene in usedGene) {
  message(gene)
  
  pdf(str_c("track/", gene, ".pdf"), width = 18, height = 5.5)
  
  grid.newpage()
  grid.draw(figList[[gene]])
  
  dev.off()
}

# promoter track plot -----------------------------------------------------

usedGene <- c("")

figList <- list()

for (gene in usedGene) {
  message(gene)
  
  gAnno <- getGeneAnno(gene)
  
  pList <- list()
  for(i in names(sel)) {
    # message(i)
    pList[[i]] <- plotTrack(
      chr = gAnno$geneChr,
      startPos = gAnno$startPos,
      endPos = gAnno$endPos,
      extendLength = 2e3,
      nbin = 2000,
      bwFile = sel[[i]][sampOrder],
      sampleName = sn,
      sampleColor = markCol[i] %>% unname %>% rep(length(sn)),
      stripColor = sampColor,
      addRegion = F,
      showName = F, showTh = F,
      ratio = 1/4)
  }
  
  pList$anno <- plotAnno(
    chr = gAnno$geneChr,
    startPos = gAnno$startPos,
    endPos = gAnno$endPos,
    extendLength = 2e3,
    gene = gene, genePos = "topleft",
    ratio = 1.5/4, dpi = 600
  )
  
  for(i in names(sel)) {
    pList[[i]]$widths[pList[[i]]$layout$l[str_detect(pList[[i]]$layout$name, "axis-l")] %>% unique] <- unit(0, "null")
    pList[[i]]$heights[pList[[i]]$layout$t[str_detect(pList[[i]]$layout$name, "axis-b")]] <- unit(0, "null")
    pList[[i]]$widths[pList[[i]]$layout$l[str_detect(pList[[i]]$layout$name, "strip-r")] %>% unique] <- unit(0.04, "null")
    pList[[i]]$widths[length(pList[[i]]$widths)] <- unit(0.02, "null")
  }
  
  pList$anno <- ggplotGrob(pList$anno)
  pList$anno$widths[pList$anno$layout$l[str_detect(pList$anno$layout$name, "axis-l")] %>% unique] <- unit(0, "null")
  pList$anno$widths[pList$anno$layout$l[str_detect(pList$anno$layout$name, "ylab-l")] %>% unique] <- unit(0, "null")
  pList$anno$heights[pList$anno$layout$t[str_detect(pList$anno$layout$name, "axis-b")] %>% unique] <- unit(0, "null")
  pList$anno$heights[pList$anno$layout$t[str_detect(pList$anno$layout$name, "xlab-b")] %>% unique] <- unit(0, "null")
  pList$anno$widths[length(pList$anno$widths)] <- unit(0.06, "null")
  
  fig <- gtable(
    widths = unit(rep(1.06, length(sel)), "null"),
    heights = unit(c(0.1 * length(sn), 0.15), "null"))
  
  for(i in 1:length(sel)) {
    fig <- gtable_add_grob(
      fig,
      grobs = pList[[i]],
      t = 1,
      l = i)
    
    fig <- gtable_add_grob(
      fig,
      grobs = pList$anno,
      t = 2,
      l = i)
  }
  
  figList[[gene]] <- fig
}

grid.newpage()
grid.draw(figList$SALL4)

for(gene in usedGene) {
  message(gene)
  
  pdf(str_c("track/", gene, "_pro.pdf"), width = 9, height = 5)
  
  grid.newpage()
  grid.draw(figList[[gene]])
  
  dev.off()
}

