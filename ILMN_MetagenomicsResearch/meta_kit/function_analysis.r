#!/usr/bin/env Rscripts
options(warn = -1)

## get work path
print_usage <- function() {
  cat('Usage:\n\tfunction_analysis.r <config_dir> <work_dir>')
  cat('Example:\n\tconfig_dir: ./config\n\twork_dir:   ./results/5.function_annotation\n')
  quit(save = 'no', status = 1)
}

get_args <- function() {
  args <- commandArgs(T)
  if (length(args) != 2) {
    print_usage()
  } else if (!dir.exists(args[1]) || !dir.exists(args[2])) {
    print_usage()  
  } else {
    return(args) 
  }
}

args <- get_args()
config.dir <- args[1]
work.dir <- args[2]

suppressMessages({
  library(tidyverse)
  library(ggsci)
  library(scales)
  library(reshape2)
  library(pheatmap)
  library(vegan)
  library(phangorn)
  library(ape)
  library(metagenomeSeq)
  library(MASS)
})

##############################  MY CONFIGURATION  ##############################
mysci.palette30 <- c(pal_npg(alpha = 1)(10), pal_npg(alpha = 0.7)(10), pal_npg(alpha = 0.3)(10))
##############################  MY CONFIGURATION  ##############################


#~ group table
print(paste("INFO", date(), paste0("function_annotation.r Work Path: ", work.dir), sep = ' - '))
group.table <- read.table(file = paste0(config.dir, '/group.txt'), header = T, sep = '\t', 
                          row.names = 1, stringsAsFactors = F, check.names = F)
group.table_variant <- group.table %>% mutate(variable = rownames(group.table))
group.count = n_distinct(group.table)
## metagenomeSeq -- prepare group file
meta.group <- data.frame(names=rownames(group.table), group=group.table$group)
meta.group_file <- paste0(config.dir, "/meta.group.txt")
write.table(meta.group, file = meta.group_file, sep='\t', row.names = F, fileEncoding = 'UTF-8')

## 1. Function Abundance 
print(paste("INFO", date(), "Function Abundance", sep = ' - '))
func_abundance_plot <- function(FILE) {
  func.table <- read.table(FILE, header = T, sep = '\t', check.names = F, stringsAsFactors = F, row.names = 1)
  func_order <- rownames(func.table)[order(rowSums(func.table))]
  func.table_order <- func.table[func_order,]
  func.table_order$func <- factor(x = func_order, levels = func_order)
  func.table_long <- melt(func.table_order, id = 'func')
  func.table_long_group <- merge(func.table_long, group.table_variant, by = 'variable')
  #~ plot
  p <- ggplot(func.table_long_group, mapping = aes(x = variable, y = value, fill = func)) +
    geom_bar(stat="identity", width = 0.8) +
    scale_fill_manual(values = mysci.palette30) +
    facet_grid( ~ group, scales = "free_x", switch = "x") +
    theme_bw() +
    labs(x = 'Sample Names',
         y = 'Relative Abundance') +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(reverse = T, ncol = 1))
  outfile <- str_replace(FILE, ".xls", ".png")
  ggsave(filename = outfile, plot = p, dpi = 300, width = 9, height = 7)
}

#~ eggNOG lv2
func.file <- paste0(work.dir, '/eggNOG/eggnog.function_lv2.merge_table.xls')
func_abundance_plot(func.file)
#~ KEGG lv1
func.file <- paste0(work.dir, '/KEGG/kegg.pathway_lv1.merge_table.xls')
func_abundance_plot(func.file)
## CAZy lv1
func.file <- paste0(work.dir, '/CAZy/cazy_func_lv1.merge_table.xls')
func_abundance_plot(func.file)

## 2. Cluster Analusis
print(paste("INFO", date(), "Cluster Analusis", sep = ' - '))
func_pheatmap_plot <- function(FILE) {
  func.table <- read.table(FILE, header = T, sep = '\t', check.names = F, stringsAsFactors = F, row.names = 1)
  func.table_head30 <- func.table %>%
    arrange(desc(rowSums(func.table))) %>%
    head(30)
  #~ draw pheatmap, add group info
  outfile <- str_replace(FILE, ".merge_table.xls", "_pheatmap.png")
  png(filename = outfile, units = "in", height = 9, width = 5, res = 300)
  pheatmap(func.table_head30, scale="column", annotation_col = group.table, cluster_cols = F)
  dev.off()
}

## 3. PCA & NMDS
func_pca_nmds_plot <- function(FILE) {
  func.table <- read.table(FILE, header = T, sep = '\t', check.names = F, stringsAsFactors = F, row.names = 1)
  func.table <- func.table %>% 
    arrange(desc(rowSums(func.table)))
  func.table <- func.table[1:50,] %>% t()
  func.rare <- as.data.frame(rrarefy(func.table, min(rowSums(func.table))))
  if (length(which(colSums(func.rare) == 0)) != 0) {func.rare <- func.rare[,-which(colSums(func.rare) == 0)]}
  func.data <- as.data.frame(func.rare)
  func.data$group <- group.table[rownames(func.data),]
  #~ PCA
  pca <- prcomp(func.data[,-ncol(func.data)], scale. = T)
  pca.var = pca$sdev^2  ## Calculate the weight of each data in the original data on each PC
  pca.var.per = round(pca.var/sum(pca.var)*100,2)  ## Calculate the ratio column of each PC to the sum of all PCs
  pca.table <- as.data.frame(pca$x) %>%
    mutate(group = func.data$group)
  p <- ggplot(aes(PC1,PC2,color = group), data = pca.table) +
    geom_point(size = 3, aes(shape=group)) +
    theme_bw() +
    labs(x = paste('PC1(',pca.var.per[1],'%)',sep = ''),
         y = paste('PC2(',pca.var.per[2],'%)',sep = ''),
         title = "PCA - PC1 vs PC2") +
    stat_ellipse(level = 0.68) + 
    theme(axis.text = element_text(size = 11,color = 'black'),
          axis.ticks = element_line(color = 'black'),
          plot.title = element_text(hjust = 0.5))
  outfile.pca <- str_replace(FILE, ".merge_table.xls", "_pca.png")
  ggsave(filename = outfile.pca, plot = p, dpi = 300, width = 7, height = 7)
  #~ NMDS
  nmds <- metaMDS(func.rare, trace = F)
  site.scores <- scores(nmds) %>% as.data.frame()
  site.scores$group <- group.table[rownames(site.scores),]
  p <- ggplot(data = site.scores,aes(x=NMDS1,y=NMDS2)) + 
    geom_point(aes(shape=group,color=group,fill=group),size=3) + 
    labs(title = "NMDS - NMDS1 vs NMDS2") + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_vline(xintercept=0,linetype="dashed",color="blue") + 
    geom_hline(yintercept = 0,linetype="dashed",color="blue") + 
    geom_text(x=0.2,y=-0.2,label=paste("stress = ",round(nmds$stress,3),sep = ""),color="blue",size=5)
  outfile.nmds <- str_replace(FILE, ".merge_table.xls", "_nmds.png")
  ggsave(filename = outfile.nmds, plot = p, dpi = 300, width = 7, height = 7)
}

## 4. PCoA, UPGMA, Anosim by Bray-Curtis distance
func_upgma_pcoa_plot <- function(FILE) {
  func.table <- read.table(FILE, header = T, sep = '\t', check.names = F, stringsAsFactors = F, row.names = 1)
  func.table <- func.table %>% arrange(desc(rowSums(func.table)))
  func.table <- func.table[1:50,] %>% t()
  func.rare <- as.data.frame(rrarefy(func.table, min(rowSums(func.table))))
  if (length(which(colSums(func.rare) == 0)) != 0) {func.rare <- func.rare[,-which(colSums(func.rare) == 0)]}
  #~ metagenomeSeq -- prepare input files
  func.meta <- as.data.frame(t(func.rare))
  Taxanomy <- rownames(func.meta)
  FUNC <- 1:nrow(func.meta)
  func.taxanomy <- data.frame(FUNC, Taxanomy)
  func.taxanomy.file <- str_replace(FILE, ".merge_table.xls", "_meta.taxa_taxanomy.csv")
  write.table(func.taxanomy, file = func.taxanomy.file, quote = T, sep='\t', row.names = F, fileEncoding = 'UTF-8')
  func.meta <- cbind(FUNC, func.meta)
  func.counts.file <- str_replace(FILE, ".merge_table.xls", "_meta.taxa_counts.csv")
  write.table(func.meta, file = func.counts.file, quote = T, sep='\t', row.names = F, fileEncoding = 'UTF-8')
  #~ Bray-Curtis distance
  distance <- vegdist(func.rare, method = 'bray')
  distance_matrix <- as.matrix(vegdist(func.rare, method = 'bray', diag=T, upper=T))
  #~ UPGMA (Bray-Curtis distance based)
  plt.png_upgma <- str_replace(FILE, ".merge_table.xls", "_upgma.png")
  up <- upgma(distance)
  png(plt.png_upgma, units = "in", height = 9, width = 9, res = 300)
  plot(up, type="phylogram", main="phylogram")
  dev.off()
  #~ PCoA
  res <- ape::pcoa(distance_matrix, correction = "none", rn = NULL)
  result <- res$values[, "Relative_eig"]
  pro1 <- as.numeric(sprintf("%.3f",result[1])) * 100
  pro2 <- as.numeric(sprintf("%.3f",result[2])) * 100
  x <- res$vectors
  sample_names <- rownames(x)
  pc <- as.data.frame(res$vectors)
  pc$names <- sample_names
  legend_title <- ""
  pc$group <- group.table[rownames(pc),]
  
  p <- ggplot(data = pc, aes(Axis.1, Axis.2, color=group, shape=group)) + 
    geom_point(size=3) + 
    stat_ellipse(level = 0.7) + 
    labs(x = paste("PCoA1(",pro1,"%)", sep = ""),
         y = paste("PCoA2(",pro2,"%)", sep = ""),
         title = "PCoA - PCoA1 vs PCoA2 (Bray Curtis)") + 
    # color=legend_title, shape=legend_title) +
    geom_hline(yintercept=0,linetype=4,color="grey") +
    geom_vline(xintercept=0,linetype=4,color="grey") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  outfile.pcoa <- str_replace(FILE, ".merge_table.xls", "_pcoa.png")
  ggsave(filename = outfile.pcoa, plot = p, dpi = 300, width = 7, height = 7)
  #~ anosim
  func.anosim <- anosim(distance_matrix, group.table[rownames(distance_matrix),], permutations=999)
  plt.png_anosim <- str_replace(FILE, ".merge_table.xls", "_anosim.png")
  png(plt.png_anosim, units = "in", height = 5, width = 5, res = 300)
  plot(func.anosim, col=mysci.palette30[1:(n_distinct(group.table)+1)], main="Anosim", xlab="", ylab="")
  dev.off()
}

## 5. metagenomeSeq
metagenomeseq_process <- function(FILE) {
  func.counts.file <- str_replace(FILE, ".merge_table.xls", "_meta.taxa_counts.csv")
  func.taxanomy.file <- str_replace(FILE, ".merge_table.xls", "_meta.taxa_taxanomy.csv")
  func.meta_counts <- read.delim(func.counts.file, row.names = 1, sep = "\t", stringsAsFactors = FALSE, 
                                 check.names = FALSE) %>% as.data.frame()
  func.meta_taxa <- data.frame(taxa=rownames(func.meta_counts))
  func.meta <- list(counts=func.meta_counts, taxa=func.meta_taxa)
  taxa <- read.delim(func.taxanomy.file, stringsAsFactors = F)
  #@ group.txt header can't startswith '#'!!!
  clin <- loadPhenoData(meta.group_file, tran = T)
  clin$siteSampled <- "intestine"  ## Placeholder, prevent data.frame --> character
  clin$group <- gsub('group', 'GROUP', clin$group)
  clin.unique <- unique(clin$group)
  rareFeatures_threshold <- length(clin$group) / 66 * 10  # origin package example threshold is '10/66' 
  
  ord <- match(colnames(func.meta$counts), rownames(clin))
  clin <- clin[ord, ]
  phenotypeData <- AnnotatedDataFrame(clin)
  OTUdata <- AnnotatedDataFrame(taxa)
  obj <- newMRexperiment(func.meta$counts, phenoData = phenotypeData, featureData = OTUdata)
  rareFeatures <- which(rowSums(MRcounts(obj) > 0) < rareFeatures_threshold)
  if ( length(rareFeatures) != 0 ) {obj <- obj[-rareFeatures, ]}  ## df[-0,]  ----  remove all rows (all data)
  obj <- filterData(obj, present = ncol(func.meta$counts)*0.375, depth = 1)
  objp <- cumNormStat(obj, pFlag = F)
  obj <- cumNorm(obj, p=objp)
  
  if (group.count == 2) {
    #~ 2 groups 
    pd <- pData(obj)
    mod <- model.matrix(~group, data = pd)
    objresl <- fitFeatureModel(obj, mod)
    metaseq.top20 <- MRcoefs(objresl, number = 20, coef = 'adjPvalues')
    metaseq.top20_simple <- metaseq.top20[,1:4]
  } else {
    #~ multiple groups 
    group <- pData(obj)$group
    settings = zigControl(maxit = 1, verbose = FALSE)
    mod = model.matrix(~-1 + group)
    colnames(mod) <- gsub("group", "", colnames(mod))
    res <- fitZig(obj, mod = mod, control = settings)
    zigFit = slot(res, "fit")
    finalMod = slot(res, "fit")$design
    contrasts <- c()
    for (a in 1:(length(clin.unique)-1)) {
      for (b in (a+1):length(clin.unique)) {
        contrast <- str_c(clin.unique[a],clin.unique[b],sep = "-")
        contrasts <- c(contrasts, contrast)}}
    contrast.matrix = makeContrasts(contrasts = contrasts, levels = finalMod)
    fit2 = contrasts.fit(zigFit, contrast.matrix)
    fit2 = eBayes(fit2)
    metaseq.top20 <- topTable(fit2, number = 20)
    metaseq.top20_simple <- metaseq.top20[,c("AveExpr","F","P.Value","adj.P.Val")]
  }
  taxids <- as.character(obj@featureData@data[rownames(metaseq.top20), ]$Taxanomy)
  metaseq.top20 <- cbind(Taxa=taxids, metaseq.top20)
  metaseq.top20_simple <- cbind(Taxa=taxids, metaseq.top20_simple)
  #~ write metagenomeseq result
  metaseq.result.file <- str_replace(FILE, ".merge_table.xls", "_metagenomeSeq.taxa_significant.xls")
  metaseq.simple_result.file <- str_replace(FILE, ".merge_table.xls", "_metagenomeSeq.taxa_significant_simple.xls")
  write.table(metaseq.top20, file = metaseq.result.file, sep = '\t', row.names = F, fileEncoding = 'UTF-8')
  write.table(metaseq.top20_simple, file = metaseq.simple_result.file, sep = '\t', row.names = F, fileEncoding = 'UTF-8')
}

## x. final plot
final_plot <- function(FILE) {
  func_pheatmap_plot(FILE)
  func_pca_nmds_plot(FILE)
  func_upgma_pcoa_plot(FILE)
  metagenomeseq_process(FILE)
}
#~ KEGG ko
func.file <- paste0(work.dir, '/KEGG/kegg.ko.merge_table.xls')
final_plot(func.file)
#~ GO go
func.file <- paste0(work.dir, '/GO/go.merge_table.xls')
final_plot(func.file)
#~ eggNOG OGs
func.file <- paste0(work.dir, '/eggNOG/eggnog.og.merge_table.xls')
final_plot(func.file)
#~ CAZy lv2
func.file <- paste0(work.dir, '/CAZy/cazy_func_lv2.merge_table.xls')
final_plot(func.file)
