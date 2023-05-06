#!/usr/bin/env Rscripts
options(warn = -1)

## get work path
print_usage <- function() {
  cat("Usage:\n\tfunction_analysis.r <config_dir> <results_dir>")
  cat("Example:\n\tconfig_dir: ./config\n\twork_dir:   ./results\n")
  quit(save = "no", status = 1)
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
results.dir <- args[2]
work.dir <- paste(results.dir, "5.function_annotation", sep = "/")

## library
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
  library(circlize)
  library(ggpubr)
})

##############################  MY CONFIGURATION  ##############################
mysci.palette30 <- c(pal_npg(alpha = 1)(10), pal_npg(alpha = 0.7)(10), pal_npg(alpha = 0.3)(10))
##############################  MY CONFIGURATION  ##############################

#~ group table
print(paste("INFO", date(), paste0("function_annotation.r Work Path: ", work.dir), sep = " - "))
group.table <- read.table(file = paste0(config.dir, "/group.txt"), header = T, sep = "\t", 
                          row.names = 1, stringsAsFactors = F, check.names = F)
group.table_variant <- group.table %>% mutate(variable = rownames(group.table))
group.count = n_distinct(group.table)
group.unique <- unique(group.table$group)

## 1. Function Abundance 
print(paste("INFO", date(), "Function Abundance", sep = " - "))
#~ aro
func.file <- paste(work.dir, "CARD/aro.merge_table.xls", sep = "/")
func.table <- read.delim(func.file, header = T, sep = "\t", check.names = F, stringsAsFactors = F, row.names = 1)
func.table <- func.table %>% 
  arrange(desc(rowSums(func.table)))
func.table_top10 <- func.table[1:10,]
func.table_top20 <- rbind(func.table[1:20,], Others=colSums(func.table[21:nrow(func.table),]))
func.table_top30 <- func.table[1:30,]
#~ aro > 3 words, translate
aro_long2short <- function(x) {
  x_vector <- strsplit(x, split = " ")[[1]]
  if (length(x_vector) > 3) {
    return(paste(x_vector[1:3], collapse = "_"))
  } else {
    return(x)
  }
}
rownames(func.table_top10) <- sapply(rownames(func.table_top10), aro_long2short)
rownames(func.table_top20) <- sapply(rownames(func.table_top20), aro_long2short)
rownames(func.table_top30) <- sapply(rownames(func.table_top30), aro_long2short)
## Top 20 ARO cumulative barplot
func.table_top20$func <- factor(x = rownames(func.table_top20), 
                                levels = rev(rownames(func.table_top20)))
func.table_top20_long <- melt(func.table_top20, id = "func")
func.table_top20_long_group <- merge(func.table_top20_long, group.table_variant, by = "variable")
#~ plot (ppm)
p <- ggplot(func.table_top20_long_group, mapping = aes(x = variable, y = value, fill = func)) +
  geom_bar(stat="identity", width = 0.8) +
  scale_fill_manual(values = mysci.palette30) +
  facet_grid( ~ group, scales = "free_x", switch = "x") +
  theme_bw() +
  labs(x = "Sample Names",
       y = "Relative Abundance(ppm)") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(reverse = T, ncol = 1))
outfile <- str_replace(func.file, ".xls", "_ppm.png")
ggsave(filename = outfile, plot = p, dpi = 300, width = 9, height = 7)
#~ plot (percent)
p <- ggplot(func.table_top20_long_group, mapping = aes(x = variable, y = value, fill = func)) +
  geom_bar(stat="identity", width = 0.8, position="fill") +
  scale_fill_manual(values = mysci.palette30) +
  facet_grid( ~ group, scales = "free_x", switch = "x") +
  theme_bw() +
  labs(x = "Sample Names",
       y = "Relative Abundance(%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = percent) +
  guides(fill = guide_legend(title = NULL, reverse = T, ncol = 1))
outfile <- str_replace(func.file, ".xls", "_percent.png")
ggsave(filename = outfile, plot = p, dpi = 300, width = 9, height = 7)

## chord diagram
print(paste("INFO", date(), "Function Abundance (chord)", sep = " - "))
color.number <- sum(dim(func.table_top10))
#~ legend map mark to ARO/sample
aro.legend <- data.frame(code=LETTERS[1:10], ARO=rownames(func.table_top10))
sample.legend <- data.frame(code=letters[1:length(colnames(func.table_top10))], sample=colnames(func.table_top10))
temp.plot_table <- func.table_top10
dimnames(temp.plot_table) <- list(aro.legend$code, sample.legend$code)
#~ need translate to matrix
temp.plot_table <- temp.plot_table %>% as.matrix()
#~ plot
outfile <- str_replace(func.file, ".xls", "_chord.png")
png(paste(outfile, sep = "/"), units="in", width=12, height=8, res=300)
par(mfrow=c(1,2))
par(mar = c(0,4,0,0) + 0.1, fig = c(0,0.6,0,1))
chordDiagram(temp.plot_table, grid.col = mysci.palette30[1:color.number], annotationTrack = c("grid", "name"))
par(fig=c(0.6,1,0,1), mai = c(0, 1, 0, 0), new =T)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,10))
legend(0,8, legend = aro.legend$ARO, pch = paste0(aro.legend$code, collapse = ""),
       col = mysci.palette30[1:10], bty="n", cex = 0.8, title = "ARO", title.adj = 0)
legend(0,5, legend = sample.legend$sample, pch =  paste0(sample.legend$code, collapse = ""),
       col = mysci.palette30[11:(10+length(sample.legend$sample))], bty="n", cex = 0.8,
       title = "Sample", title.adj = 0)
circos.clear()
dev.off()
## anosim
print(paste("INFO", date(), "Anosim Plot", sep = " - "))
func.table_t <- t(func.table)
func.anosim <- anosim(func.table_t, group.table[rownames(func.table_t),], permutations=999)
outfile <- str_replace(func.file, ".xls", "_anosim.png")
png(outfile, units = "in", height = 5, width = 5, res = 300)
plot(func.anosim, col=mysci.palette30[1:(n_distinct(group.table)+1)], main="Anosim", xlab="", ylab="")
dev.off()

## pheatmap
print(paste("INFO", date(), "Pheatmap", sep = " - "))
#~ draw pheatmap, add group info
outfile <- str_replace(func.file, ".xls", "_pheatmap.png")
png(filename = outfile, units = "in", height = 9, width = 6, res = 300)
pheatmap(func.table_top30, scale="column", annotation_col = group.table, cluster_cols = F)
dev.off()

## difference - boxplot & pointplot
fig.table <- cbind(group.table, aro_count = colSums(func.table > 0)[rownames(group.table)],
                   sample=rownames(group.table))
#~ pairwise comparison between groups
if (group.count == 2) {
  comparison <- list(group.unique)
} else {
  comparison <- list()
  for (f in 1:(group.count-1)) {
    for (b in (f+1):group.count) {comparison <- append(comparison,list(c(group.unique[f], group.unique[b])))}
  }
}
#~ label y coeficient
if (group.count == 2) {
  label.y_coef <- 1.1
} else {
  label.y_coef <- 1.5
}
#~ label y value
rng <- range(fig.table["aro_count"], na.rm = T)
label.y_value <- rng[1] + diff(rng) * label.y_coef

p <- ggboxplot(fig.table, x = "group", y = "aro_count", color = "group",
               palette = mysci.palette30, add = "jitter", shape = "group") +
  stat_compare_means(label.y = label.y_value) +
  stat_compare_means(comparisons = comparison, label = "p.signif") +
  labs(x = "Groups", y = "Number of AROs")
outfile <- outfile <- str_replace(func.file, ".xls", "_boxplot.png")
ggsave(filename = outfile, plot = p, dpi = 300)

## doughnut chart : aro species composition VS kraken species composition
for (grp in group.unique) {
  temp_group.table <- filter(group.table, group == grp)
  #~ multiple groups
  aro_species_whole.table <- matrix(ncol = 1) %>% as.data.frame()
  colnames(aro_species_whole.table) <- c("phylum")
  kraken_annotation.table <- aro_species_whole.table
  for (samp in rownames(temp_group.table)) {
    aro_igi_species.file <- paste(work.dir, samp, "aro_rgi_simple.xls", sep = "/")
    kraken_annotation.file <- paste(results.dir, "4.species_annotation", samp, 
                                    "gene_annotation_by_kraken2.txt", sep= "/")
    aro_temp.table <- read.delim(aro_igi_species.file, header = T, sep = "\t", check.names = F,
                                 stringsAsFactors = F)['phylum']
    kraken_temp.table <- read.delim(kraken_annotation.file, header = T, sep = "\t", check.names = F,
                                    stringsAsFactors = F)['phylum']
    aro_species_whole.table <- rbind(aro_species_whole.table, aro_temp.table)
    kraken_annotation.table <- rbind(kraken_annotation.table, kraken_temp.table)
  }
  aro_species_whole_count.table <- aro_species_whole.table %>% 
    filter(!is.na(phylum) & phylum != "None") %>% 
    group_by(phylum) %>% 
    summarise(aro_count = n())
  kraken_species_whole_count.table <- kraken_annotation.table %>% 
    filter(!is.na(phylum) & phylum != "None") %>% 
    group_by(phylum) %>% 
    summarise(kraken_count = n())
  #~ merge kraken aro table
  group_merge.table <- merge(kraken_species_whole_count.table, aro_species_whole_count.table, 
                             by = "phylum", all = T)
  group_merge.table[is.na(group_merge.table)] <- 0
  group_merge.table <- group_merge.table %>% arrange(desc(rowSums(group_merge.table[,-1])))
  group_merge.table <- rbind(group_merge.table[1:3,],
        c("Others", colSums(group_merge.table[4:dim(group_merge.table)[1],-1]))) %>% 
    mutate_at(c("kraken_count", "aro_count"), as.numeric) %>% 
    #~ label
    mutate(kraken_percent = paste0(round(kraken_count / sum(kraken_count), digits = 2) * 100, "%")) %>% 
    mutate(aro_percent = paste0(round(aro_count / sum(aro_count), digits = 2) * 100, "%")) %>% 
    #~ fraction
    mutate(kraken_fraction = kraken_count / sum(kraken_count)) %>%
    mutate(aro_fraction = aro_count / sum(aro_count)) %>% 
    #~ cumulative percentages, top of each rectangle
    mutate(kraken_ymax = cumsum(kraken_fraction)) %>% 
    mutate(aro_ymax = cumsum(aro_fraction)) %>% 
    #~ bottom of each rectangle
    mutate(kraken_ymin = c(0, head(cumsum(kraken_fraction), n=-1))) %>% 
    mutate(aro_ymin = c(0, head(cumsum(aro_fraction), n=-1))) %>% 
    #~ label position
    mutate(kraken_labelPosition = (kraken_ymax + kraken_ymin) / 2) %>% 
    mutate(aro_labelPosition = (aro_ymax + aro_ymin) / 2) 
  group_merge.table$phylum <- factor(group_merge.table$phylum, levels = rev(group_merge.table$phylum))
  #~ inner : ARO,  outer : kraken(whole)
  p <- ggplot(group_merge.table) +
    geom_col(aes(x = 2, y = aro_count, fill = phylum), position = "fill", col = "white") +
    geom_col(aes(x = 3, y = kraken_count, fill = phylum), position = "fill", col = "white") +
    geom_text(aes(x = 2, y = aro_labelPosition, label = aro_percent)) +
    geom_text(aes(x = 3, y = kraken_labelPosition, label =kraken_percent)) +
    xlim(0, 3.5) + 
    labs(x = NULL, y = NULL, title = grp) + 
    coord_polar(theta = "y") +
    scale_fill_manual(values = mysci.palette30) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank()) +
    guides(fill = guide_legend(title = NULL, reverse = T))
  out.file <- str_replace(func.file, ".xls", paste0("_", grp, "_doughnut.png"))
  ggsave(filename = out.file, plot = p, dpi = 600)
}
