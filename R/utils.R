
setCMD <- function(cmd, dir = ".", sepN = 1) {
  cmd %>% tapply(seq_along(.) %% sepN, c) %>% map(~ {
    c("#!/bin/bash", .x)}) %T>%
    iwalk(~ write.table(.x, glue("{dir}/batch{.y}.sh"), quote = F, row.names = F, col.names = F)) %>%
    names() %>% map_chr(~ glue(
      "{head} {dir}/batch{.x}.sh {tail}",
      head = "sh", tail = "&")) %>%
    c("#!/bin/bash", .) %>% write.table(glue("{dir}/submit.sh"), quote = F, row.names = F, col.names = F)
}

plotATACmeta <- function(
    obj, fs = "Sample", fs_type = c("d", "c"), lb_fs = fs, reduc = "UMAP", dims = c(1, 2),
    col = pal_d3("category20")(20), lb_col = col,
    c_col = brewer.pal(11, "Spectral") %>% rev,
    x_name = str_c(toupper(reduc), "_", dims[1]), y_name = str_c(toupper(reduc), "_", dims[2]),
    size = .5, alpha = 1, dpi = 300, order_by = NULL,
    order = F, label = T
){
  require(data.table)
  require(tidyverse)
  require(ggplot2)
  require(ggrastr)
  require(RColorBrewer)
  require(ggsci)
  
  pData <- getEmbedding(obj, reduc) %>% as.data.table(keep.rownames = T) %>% set_colnames(c("rn", "dim1", "dim2"))
  pData[, meta := getCellColData(obj, fs)[[1]]]
  pData[, lb := getCellColData(obj, lb_fs)[[1]]]
  # m_x <- min(pData$dim1) - .05*diff(range(pData$dim1))
  # m_y <- min(pData$dim2) - .05*diff(range(pData$dim2))
  # q_x <- m_x + .2*diff(range(pData$dim1))
  # q_y <- m_y + .2*diff(range(pData$dim2))
  
  if(fs_type[1] == "c") {
    colorBar <- scale_color_gradientn(colors = c_col)
    
    guideUse <- guides()
    
    if(order) {pData <- pData[order(meta)]} else {pData <- pData[sample(1:.N)]}
  } else if(fs_type[1] == "d") {
    n <- length(unique(pData$meta))
    if(length(col) < n) {col <- colorRampPalette(col)(n)}
    colorBar <- scale_color_manual(values = col)
    
    guideUse <- guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 1)))
    
    if(order) {pData <- pData[order_by]} else {pData <- pData[sample(1:.N)]}
  }
  
  p <- ggplot(pData, aes(x = dim1, y = dim2)) +
    geom_point_rast(aes(color = meta), size = size, alpha = alpha, shape = 16, raster.dpi = dpi) +
    # geom_segment(x = m_x, y = m_y, xend = q_x, yend = m_y) +
    # geom_segment(x = m_x, y = m_y, xend = m_x, yend = q_y) +
    colorBar +
    guideUse +
    labs(x = x_name, y = y_name, title = fs, color = "") +
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray30"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(hjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(face = "bold", hjust = .5))
  
  if(fs_type == "d" && label){
    label_x <- tapply(pData$dim1, pData$lb, median)
    label_y <- tapply(pData$dim2, pData$lb, median)
    label_text <- names(label_x)
    
    lData <- data.frame(
      x = label_x,
      y = label_y,
      label = label_text)
    p <- p +
      geom_label(
        data = lData, aes(x = x, y = y, label = label),
        fill = "white", size = 3, color = NA, alpha = .6, fontface = "bold",
        label.size = 0.5, label.r = unit(0.25, "lines"), label.padding = unit(0.15, "lines")) +
      geom_text(
        data = lData, aes(x = x, y = y, label = label),
        size = 3, color = "black", fontface = "bold")
  }
  p
}

getMatColor <- function(pal = "blue", num = 3) {
  pal_material(palette = pal, n = 2 * num + 1)(2 * num + 1)[2 * 1:num]
}

plotRectVenn <- function(
    pct = 0.5, txt = c("a", "b", "over"),
    label = c("group1", "group2"),
    typeCol = pal_d3()(3)
) {
  fei <- (5^.5 - 1)/2
  # dx <- (fei+1 - ((fei+1)^2 - 4*fei*(1-pct))^.5)/(2*fei)
  dx <- 1 - pct
  
  (1 - dx)*(1 - fei*dx)
  
  rectData <- data.table(
    xmin = c(0, dx, dx),
    xmax = c(1, 1 + dx, 1),
    ymin = c(0, fei * dx, fei * dx),
    ymax = c(1, 1 + fei * dx, 1),
    type = c("a", "b", "over")
  )
  
  textData <- data.table(
    x = c(dx/2,  1 + dx/2, 1/2 + dx/2),
    y = c(1/2, 1/2 + fei * dx, 1/2 + 1/2 * fei * dx),
    label = txt
  )
  
  labelData <- data.table(
    x = c(1/2, 1/2 + dx),
    y = c(-0.15, 1.15 + fei * dx),
    type = c("a", "b"),
    label = label
  )
  
  ggplot(rectData, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)) +
    geom_rect(aes(fill = type), alpha = 1, color = "black", show.legend = F) +
    geom_text(data = textData, aes(x = x, y = y, label = label), size = 5) +
    geom_text(
      data = labelData, aes(x = x, y = y, label = label, color = type),
      size = 6, fontface = "bold", show.legend = F) +
    scale_fill_manual(values = typeCol) +
    scale_color_manual(values = typeCol) +
    scale_x_continuous(limits = c(0, 2)) +
    scale_y_continuous(limits = c(-0.3, 2.3), expand = expansion(c(0, 0))) +
    coord_fixed() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}

peakL <- function(x) {
  str_split(x, "_") %>% map_dbl(~ as.numeric(.x[3]) - as.numeric(.x[2]))
}

plotTrack <- function(
    gene = "SALL4", extendLength = 2000, nbin = 1e3, binSize = 20, useBS = F,
    ratio = 1/10, dpi = 300,
    hasSig = T, minTh = 0, sizeFac = rep(1, length(bwFile)), bgSig = rep(0, length(bwFile)),
    sampColor = ArchR::ArchRPalettes$stallion %>% unname,
    geneInfo = getGeneAnno(gene),
    addRegion = F, selRegion,
    bwFile, addChr = F, sampleName = str_remove(bwFile, ".bw$"),
    showName = T, showTh = T
) {
  geneChr <- geneInfo$geneChr
  if(addChr) {geneChr <- str_c("chr", geneChr)}
  startPos <- geneInfo$startPos
  endPos <- geneInfo$endPos
  
  extend <- max(extendLength, (endPos - startPos) * 0.05)
  startEx <- startPos - extend
  endEx <- endPos + extend
  
  usedData <- map(bwFile, ~ import.bw(
    .x, format = "bw", selection = BigWigSelection(GRanges(geneChr, IRanges(startEx, endEx)))))
  names(usedData) <- sampleName
  
  usedData %<>% map(~ as.data.table(.x))
  
  usedData %<>% imap(~ {
    posData <- data.table(
      pos = .x[1, start]:.x[.N, end],
      value = rep(.x$score, .x$width)
    )
    if(useBS) {posData[, bin := round(pos / binSize) * binSize]} else {
      posData[, bin := cut_interval(pos, nbin) %>% as.numeric]
    }
    binData <- posData[, .(sample = .y, value = mean(value)), by = bin]
    binData %>% set_colnames(c("pos", "sample", "value"))
  })
  
  # for geom_area
  # usedData %<>% imap(~ {
  #   posData <- data.table(pos = .x[, c("start", "end")] %>% as.matrix %>% t %>% as.vector, value = rep(.x$score, each = 2))
  #   if(useBS) {posData[, bin := round(pos / binSize) * binSize]} else {
  #     posData[, bin := cut_interval(pos, nbin) %>% as.numeric]
  #   }
  #   binData <- posData[, .(sample = .y, value = mean(value)), by = bin]
  #   binData %>% set_colnames(c("pos", "sample", "value"))})
  
  if(!useBS) {binSize <- (endEx - startEx)/nbin}
  bgSig <- bgSig * binSize / 20
  usedData <- pmap(list(usedData, bgSig, sizeFac), function(x, y, z) {x[, value := (value - y)/z][value < 0, value := 0]})
  
  plotData <- do.call(rbind, usedData)
  plotData$sample %<>% factor(levels = sampleName)
  
  valMax <- max(plotData$value)
  if(!hasSig && max(plotData$value) < minTh) {
    valMax <- minTh
  }
  
  if(addRegion) {
    reg <- geom_rect(
      data = selRegion,
      mapping = aes(xmin = start, xmax = end), ymin = 0, ymax = valMax * 1.05,
      fill = rep(selRegion$color, length(sampleName)), color = NA)
  } else {
    reg <- NULL
  }
  
  usedTheme <- theme(
    aspect.ratio = ratio,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    panel.spacing = unit(0, "cm"),
    panel.spacing.y = unit(0, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    plot.margin = margin(),
    axis.ticks.y = element_line(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(),
    axis.title = element_blank()
  )
  
  if(!showName) {
    usedTheme$strip.text <- element_blank()
  }
  
  if(!showTh) {
    usedTheme$axis.ticks.y <- element_blank()
    usedTheme$axis.text.y <- element_blank()
  }
  
  ggplot(plotData) +
    reg +
    rasterise(geom_col(aes(fill = sample, x = pos, y = value), width = 1, show.legend = F), dpi = dpi) +
    scale_fill_manual(values = sampColor) +
    labs(y = "Signal intensity") +
    scale_y_continuous(
      expand = c(0, 0), limits = c(0, valMax * 1.05),
      breaks = c(valMax), labels = round(valMax, 2)) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_wrap( ~ sample, ncol = 1, strip.position = "right") +
    usedTheme
}

plotAnno <- function(
    gene = "SALL4", extendLength = 2000,
    geneInfo = getGeneAnno(gene),
    ratio = 1/10
) {
  startPos <- geneInfo$startPos
  endPos <- geneInfo$endPos
  strand <- geneInfo$strand
  exonData <- geneInfo$exonData
  
  segmentData <- data.table(
    start = seq(startPos, endPos, length.out = 10)[2:9]
  )
  
  if(strand == "-") {
    segmentData[, end := start - 1]
  }else{
    segmentData[, end := start + 1]
  }
  
  geneRange <- endPos - startPos
  sizeBar <- log10(geneRange) %>% floor() %>% 10^.
  
  extend <- max(extendLength, geneRange * 0.05)
  startEx <- startPos - extend
  endEx <- endPos + extend
  
  ggplot(exonData) +
    geom_linerange(x = 0, ymin = startPos, ymax = endPos, size = 0.5) +
    geom_linerange(aes(x = 0, ymin = start, ymax = end), size = 3) +
    geom_segment(
      data = segmentData, aes(y = start, yend = end), x = 0, xend = 0,
      arrow = arrow(length = unit(.4, "lines")), size = 0.3) +
    annotate("segment", x = 1, xend = 1, y = endEx - sizeBar, yend = endEx, size = 1.5) +
    annotate("text", x = 1.2, y = endEx, label = humanFormat(sizeBar), hjust = 1, size = 3) +
    scale_y_continuous(expand = c(0, 0), limits = c(startEx, endEx)) +
    scale_x_continuous(expand = expansion(c(.1, .1))) +
    labs(x = "", y = gene) +
    coord_flip() +
    theme(
      aspect.ratio = ratio,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(face = "bold.italic"),
      plot.margin = margin()
    )
}

humanFormat <- function(n) {
  stopifnot(is.numeric(n), n > 0)
  
  size <- log10(n)
  
  if(size < 3) {
    as.character(n)
  } else if(size < 6) {
    paste0(round(n / 1e3), "K")
  } else if(size < 9) {
    paste0(round(n / 1e6), "M")
  } else {
    paste0(round(n / 1e9), "G")
  }
}

