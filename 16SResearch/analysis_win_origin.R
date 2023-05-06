options(warn = -1)

suppressMessages({
  library(parallel)
  library(vegan)
  library(tidyverse)
  library(pheatmap)
  library(scales)
  library(reshape2)
  library(RColorBrewer)
  library(phangorn)
  library(FSA)  ## dunnTest
  library(ggalt)
  library(ggsci)
  library(stringi)
  library(ggpubr)
  library(xlsx)
  library(rjson)
  library(metagenomeSeq)
  library(VennDiagram)
  })

###########################  MY CONFIGURATION  ################################
my.wilcoxon <- function(x, d=otu.use) {
  wt <- wilcox.test(x ~ d$group, data = d, exact = F)
  return(c(wt$statistic, wt$p.value)) 
}

my.kruskal <- function(x, d=otu.use) {
  kw <- kruskal.test(x ~ d$group, data = d)
  return(c(kw$statistic, kw$p.value))
}

my.dunntest <- function(x, d=alpha.table) {
  d$group <- as.factor(d$group)
  dunnTest(x ~ d$group, data = d, method="bonferroni")
}

my.ttest <- function(x, d=otu.use) {
  tt <- t.test(x ~ d$group, data = d)
  tt.temp <- c(tt$statistic, tt$p.value)
  names(tt.temp) <- c("t", "p.value")
  return(tt.temp)
}

my.aov <- function(x, d=otu.use) {
  aov.summary <- summary(aov(x ~ d$group, data = d))
  av <- aov.summary[[1]][1,4:5] %>% as.numeric()
  names(av) <- c("F value", "Pr(>F)")
  return(av)
}

mysci.palette30 <- c(pal_npg(alpha = 1)(10), pal_npg(alpha = 0.7)(10), pal_npg(alpha = 0.3)(10))

## plot output pdf & png
plot.pdf_png0 <- function(plt.pdf, plt.png, p) {
  pdf(plt.pdf)
  eval(parse(text = p))
  dev.off()
  png(plt.png)
  eval(parse(text = p))
  dev.off()
}
plot.pdf_png <- function(plt.pdf, plt.png, p, pdf_h = 9, pdf_w = 9, png_h = 720, png_w = 720) {
  pdf(plt.pdf, height = pdf_h, width = pdf_w)
  eval(parse(text = p))
  dev.off()
  png(plt.png, height = png_h, width = png_w)
  eval(parse(text = p))
  dev.off()
}
###############################################################################

## result path
result.path <- 'D:/mengxf/2021/01/LYL/data/result/sample_group_1'
OTU_analysis.path <- str_c(result.path,"/1.Taxa_analysis/")
Alpha_diversity.path <- str_c(result.path,"/2.Alpha_diversity/")
Beta_diversity.path <- str_c(result.path,"/3.Beta_diversity/")
Intermediate.path <- str_c(result.path,"/x.Intermediate/")
## confirm directory exists
if (!dir.exists(OTU_analysis.path)) {
dir.create(OTU_analysis.path, recursive = T)
dir.create(Alpha_diversity.path, recursive = T)
dir.create(Beta_diversity.path, recursive = T)
dir.create(Intermediate.path, recursive = T)
}

## otu table
raw.data <- read.xlsx2('./data/merge.xls', 1, stringsAsFactors = F, row.names = 1)
otu.table1 <- apply(raw.data, 2, as.numeric) %>% t() %>% as.data.frame()
colnames(otu.table1) <- rownames(raw.data)
## reads count threshold
otu.table2 <- otu.table1[-c(which(rowSums(otu.table1) <= 301)),]
otu.table3 <- otu.table2[,-c(which(colSums(otu.table2) == 0))]
## group 
otu.group <- read.table('./分组/groups_dir/sample_group_1.txt', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
group.unique <- unique(otu.group$group)
group.count <- length(group.unique)
## otu 
otu.table3$group <- otu.group[rownames(otu.table3),]
otu.table4 <- filter(otu.table3, !is.na(group))
otu <- otu.table4[,-ncol(otu.table4)]
sample.names <- rownames(otu)
sample.count <- length(sample.names)
## metagenomeSeq group file 
meta.group <- data.frame(names=rownames(otu.table4), group=otu.table4$group)
meta.group_file <- paste0(Intermediate.path, "meta.group.txt")
write.table(meta.group, file = meta.group_file, sep='\t', row.names = F, fileEncoding = 'UTF-8')

## real step
min.depth <- min(colSums(otu))
if (min.depth < 1000) {real_step <- 100 } else {real_step <- signif(min.depth/6, digits = 2)}

## rarefy first ## 
otu.rare <- as.data.frame(rrarefy(otu, min(rowSums(otu))))
allzero.cols <- which(colSums(otu.rare) == 0) %>% as.numeric()
if (length(allzero.cols) != 0) {otu.rare <- otu.rare[, -allzero.cols]}
rarefy.file <- paste0(OTU_analysis.path, "taxa_rarefy.txt")
write.table(t(otu.rare), file=rarefy.file, quote = T, sep="\t", row.names = T, col.names = T, fileEncoding = "UTF-8")
## metagenomeSeq -- prepare input files
otu.meta <- as.data.frame(t(otu.rare))
Taxanomy <- rownames(otu.meta)
OTU <- 1:nrow(otu.meta)
otu.taxanomy <- data.frame(OTU, Taxanomy)
otu.taxanomy.file <- str_c(Intermediate.path, "meta.taxa_taxanomy.csv")
write.table(otu.taxanomy, file = otu.taxanomy.file, quote = T, sep='\t', row.names = F, fileEncoding = 'UTF-8')
otu.meta <- cbind(OTU, otu.meta)
otu.counts.file <- str_c(Intermediate.path, "meta.taxa_counts.csv")
write.table(otu.meta, file = otu.counts.file, quote = T, sep='\t', row.names = F, fileEncoding = 'UTF-8')

#### Intra Group (alpha diversities) 
## compute Alpha Diversities Indices, result returned to vector
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)	#richness index
  else if (method == 'chao1') result <- estimateR(x)[2, ]	#Chao1 index
  else if (method == 'ace') result <- estimateR(x)[4, ]	#ACE index
  else if (method == 'shannon') result <- vegan::diversity(x, index = 'shannon', base = base)	#Shannon index
  else if (method == 'simpson') result <- vegan::diversity(x, index = 'simpson')	#Gini-Simpson index
  else if (method == 'invsimpson') result <- vegan::diversity(x, index = 'invsimpson')	#Inverse Gini-Simpson index
  else if (method == 'pielou') result <- vegan::diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou uniformity
  else if (method == 'goods_coverage') result <- 1 - rowSums(x == 1) / rowSums(x)	#goods_coverage
  else if (method == 'pd' & !is.null(tree)) {	#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

## According to the sampling step, the Alpha diversity index under each dilution gradient is calculated, and the results are returned to the list
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}

abundance <- apply(otu.rare, 1, function(x) {x/sum(x)})
richness <- alpha_index(otu.rare, method = 'richness')
goods_coverage <- alpha_index(otu.rare, method = 'goods_coverage')
chao1 <- alpha_index(otu.rare, method = 'chao1')
ace <- alpha_index(otu.rare, method = 'ace')
shannon <- alpha_index(otu.rare, method = 'shannon')
simpson <- alpha_index(otu.rare, method = 'simpson')
invsimpson <- alpha_index(otu.rare, method = 'invsimpson')
pielou <- alpha_index(otu.rare, method = 'pielou')
alpha.table <- rbind(richness, goods_coverage, chao1, ace, simpson, invsimpson, shannon, pielou) %>% t() %>% as.data.frame()
alpha.table$group <- otu.group[rownames(alpha.table),]
alpha.tablecb <- cbind(sample=rownames(alpha.table), alpha.table)
alpha.file <- paste0(Alpha_diversity.path, "alpha_table.xls")
write.table(alpha.tablecb, file=alpha.file, quote = F, sep="\t",row.names = F, fileEncoding = "UTF-8")

## alpha indices boxplot 
alpha.index_list <- c("richness", "chao1", "ace", "shannon", "simpson", "invsimpson", "pielou", "goods_coverage")
comparison <- list()
if (group.count == 2) {comparison <- list(group.unique)} else {
  comparison <- list()
  for (f in 1:(group.count-1)) {
    for (b in (f+1):group.count) {comparison <- append(comparison,list(c(group.unique[f], group.unique[b])))}
  }
}
if (group.count == 2) {label.y_coef <- 1.1} else {label.y_coef <- 1.5}
# if (group.count == 2) {ggboxplot_method <- 't.test'} else {ggboxplot_method <- 'anova'}
#~ plot 
for (ai in alpha.index_list) {
  rng <- range(alpha.table[ai], na.rm = T)
  label.y_value <- rng[1] + (rng[2] - rng[1]) * label.y_coef
  p <- ggboxplot(alpha.table, x = "group", y = ai, color = "group",
               palette = mysci.palette30, add = "jitter", shape = "group") +
  stat_compare_means(label.y = label.y_value) +
  stat_compare_means(comparisons = comparison, label = "p.signif")
ggsave(str_c(ai,"_boxplot.png"), path = Alpha_diversity.path, plot = p, dpi = 300)
}

## alpha indices difference
my.kruskal2 <- function(x) my.kruskal(x, d=alpha.table)
if ( group.count < 3 ) {
  ## K-W test
  alpha.diff <- apply(alpha.table[,-ncol(alpha.table)], 2, my.kruskal2)
  alpha.diff <- rbind(alpha.diff, p.adjust(alpha.diff[2,], method = "fdr"))
  rownames(alpha.diff)[2:3] <- c("p.value", "p.adjust")
  # alpha.diff <- cbind(Stat=rownames(alpha.diff), alpha.diff)
  alpha.diff <- cbind(Stat=colnames(alpha.diff), t(alpha.diff))
  # print(alpha.diff)
  diff.file <- paste0(Alpha_diversity.path, "alpha_diff.xls")
  write.table(alpha.diff, file = diff.file, quote = F, sep='\t', row.names = F, col.names = T, fileEncoding = "UTF-8")
} else {
  dunns <- apply(alpha.table[,-ncol(alpha.table)], 2, my.dunntest)
  for (di in 1:length(dunns)) {
    dunn.index <- names(dunns)[di]
    dunn.temp <- as.data.frame(dunns[[di]]$res)
    diff.file <- paste0(Alpha_diversity.path, dunn.index, "_alpha_diff.xls")
    write.table(dunn.temp, file = diff.file, quote = F, sep="\t",row.names = F, col.names = T, fileEncoding = "UTF-8")
  }
}

## Alpha indices dilution curve
plot_alpha <- function(method = "richness") {
  alpha.curves <- alpha_curves(otu.rare, step = real_step, method = method)
  #get ggplot2 picture file
  plot.alpha <- data.frame()
  for (i in names(alpha.curves)) {
    alpha.curves_i <- (alpha.curves[[i]])
    alpha.curves_i <- data.frame(rare = names(alpha.curves_i), alpha = alpha.curves_i, sample = i, stringsAsFactors = FALSE)
    plot.alpha <- rbind(plot.alpha, alpha.curves_i)
  }
  rownames(plot.alpha) <- NULL
  plot.alpha$rare <- as.numeric(plot.alpha$rare)
  plot.alpha$alpha <- as.numeric(plot.alpha$alpha)
  #ggplot2 plot
  p <- ggplot(plot.alpha, aes(rare, alpha, color = sample)) +
    geom_line() + geom_point() + 
    labs(x = 'Number of sequences', y = stri_trans_totitle(method), color = NULL) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(0, max(rowSums(otu.rare)), real_step), labels = as.character(seq(0, max(rowSums(otu.rare)), real_step)))
  ggsave(str_c(stri_trans_totitle(method),"_Curves.png"), path = Alpha_diversity.path, plot = p, height = 6, width = 12)
}

block <- sapply(alpha.index_list, plot_alpha)

## Top 20 species abundance
## Build taxonomy scientific map function
otu.use = t(otu.rare)
if(dim(otu.use)[1] > 20){
  top20_taxid = names(tail(sort(rowSums(otu.use), decreasing = F),20))
  Other = colSums(otu.use[!(rownames(otu.use) %in% top20_taxid),])
  otu.bar = rbind(Other, otu.use[top20_taxid,])
  otu.bar = t(t(otu.bar)/colSums(otu.bar)) * 100
}else{
  top20_taxid = names(sort(rowSums(otu.use), decreasing = F),20)
  otu.bar = otu.use[top20_taxid,]
  otu.bar = t(t(otu.use)/colSums(otu.use)) * 100
}
otu.bar <- data.frame(otu.bar, check.names = F)
# otu.bar$latin <- as.character(sapply(rownames(otu.bar), tn.func))
otu.bar$latin <- rownames(otu.bar)
otu.bar$latin <- factor(otu.bar$latin, levels = otu.bar$latin)

otu.top_long <- melt(otu.bar, id="latin")
group.temp <- data.frame(variable=sample.names, group=otu.group[sample.names,])
otu.top_long <- merge(otu.top_long, group.temp, by = 'variable')

p <- ggplot(otu.top_long, aes(x=variable, y=value, fill=latin)) + 
  geom_bar(position="fill", stat="identity", width = 0.8) +
  scale_y_continuous(labels = percent) +  ## Percentage coordinate axis
  scale_fill_manual(values = mysci.palette30) +
  theme_bw() +
  facet_grid( ~ group, scales = "free_x", switch = "x") +
  labs(x="Sample",y="Relative Abundance", title="Top 20 SPecies Abundance") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(family = 'SimSun'),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right", legend.title = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE))
# guides(fill = guide_legend(reverse = TRUE, ncol = 11, byrow = T))
ggsave("Top_20_SPecies_Abundance.png", path = OTU_analysis.path, width = 18, height = 9)

## multi-sample Rank-abundance
abundance.m <- data.frame()
for (n in 1:sample.count) {
  abundance.sub <- data.frame(sort(abundance[,n], decreasing = TRUE))
  colnames(abundance.sub) <- 'value'
  abundance.sub <- subset(abundance.sub, abundance.sub$value > 0)
  sub_rank <- 1:lengths(abundance.sub)
  abundance.sub$rank <- sub_rank
  abundance.sub$sample <- colnames(abundance)[n]
  abundance.m <- rbind(abundance.m, abundance.sub)
}

p <- ggplot(data = abundance.m, aes(x=rank, y=value, group=sample, color=sample)) + geom_line() + 
  labs(x="Taxa Rank", y="Relative Abundance", title="Rank-abundance") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(labels = scientific, limits = c(1e-3, 1e-1), trans = "log10")
ggsave("Rank_abundance.png", path = Alpha_diversity.path)

## Species Accumulation Curves
all <- specaccum(otu, method = "random")
plt.pdf <- paste0(Alpha_diversity.path, "Species_Accumulation_Curves.pdf")
plt.png <- paste0(Alpha_diversity.path, "Species_Accumulation_Curves.png")
p <- 'plot(all, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue", main = "Species Accumulation Curves", xlab = "Number of samples", ylab = "Number of  species")
boxplot(all, col = "yellow", add = TRUE, pch = "+")'
plot.pdf_png0(plt.pdf = plt.pdf, plt.png = plt.png, p = p)

## taxa distribution venn
venn.otu_rare <- otu.rare %>% mutate(group=otu.group[rownames(otu),])
venn.otu_gather <- gather(venn.otu_rare, key = "taxa", value = "count", -group) %>% filter(count>0) %>% group_by(group, taxa) %>% summarise(count=sum(count, na.rm = T))

venn.list <- list()
for (gp in group.unique) {
  venn.list[gp] <- list(filter(venn.otu_gather, group==gp)$taxa)
}
venn.file <- paste0(OTU_analysis.path, "venn.png")
venn.diagram(venn.list, filename = venn.file, 
             fill = pal_npg()(length(venn.list)), cex = 1, alpha = 0.5,
             col = 'transparent', fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif', margin = 0.2)


#### Inter (beta diversities) 
## difference analysis between groups 
## compute Bray-Curtis distance & distance matrix heatmap
distance <- vegdist(otu.rare, method = 'bray')
distance_matrix <- as.matrix(vegdist(otu.rare, method = 'bray', diag=T, upper=T))
rownames(distance_matrix) = colnames(distance_matrix) = rownames(otu.rare)
dist_matrix.file <- str_c(Beta_diversity.path, "Bary_distance.txt")
write.table(distance_matrix, file = dist_matrix.file, quote = F, sep="\t",row.names = T, col.names = T, fileEncoding = "UTF-8")

## BC distance pheatmap
plt.pdf <- paste0(Beta_diversity.path, "Bray_Curtis_distance_pheatmap.pdf")
plt.png <- paste0(Beta_diversity.path, "Bray_Curtis_distance_pheatmap.png")
p <- 'pheatmap(distance_matrix, main = "Bray-Curtis distance")'
plot.pdf_png(plt.pdf = plt.pdf, plt.png = plt.png, p = p)

## UPGMA (Bray-Curtis distance based)  # library(phangorn)
plt.pdf <- paste0(Beta_diversity.path, "upgma.pdf")
plt.png <- paste0(Beta_diversity.path, "upgma.png")
up <- upgma(distance)
p = 'plot(up, type="phylogram", main="phylogram")'
plot.pdf_png(plt.pdf = plt.pdf, plt.png = plt.png, p = p, png_h = 900, png_w = 900, pdf_h = 11, pdf_w = 11)

## PCA & adonis(PERMANOVA)
otu.data <- as.data.frame(otu.rare)
otu.data$group <- otu.group[rownames(otu.data),]
pca <- prcomp(otu.data[,-ncol(otu.data)], scale. = T)
pca.var = pca$sdev^2  ## Calculate the weight of each data in the original data on each PC
pca.var.per = round(pca.var/sum(pca.var)*100,2)  ## Calculate the ratio column of each PC to the sum of all PCs

site <- data.frame(sample = rownames(otu.data), group = otu.data$group)
otu.adonis <- adonis2(distance ~ group, site, permutations = 999)

pca.table <- as.data.frame(pca$x) %>%
  mutate(group = otu.data$group)

adonis_label <- pca.table %>% 
  summarise(
    PC1 = max(PC1),
    PC2 = max(PC2),
    group = max(group),
    label = paste("ADONIS", "Bray-Curtis", str_c("R2 = ",round(otu.adonis$R2[1],3)), str_c("Pr(>F) = ",otu.adonis$`Pr(>F)`[1]), sep = "\n"))

p <- ggplot(aes(PC1,PC2,color = group), data = pca.table) +
  geom_point(size = 3, aes(shape=group)) +
  theme_bw() +
  labs(x = paste('PC1(',pca.var.per[1],'%)',sep = ''),
       y = paste('PC2(',pca.var.per[2],'%)',sep = ''),
       title = "PCA - PC1 vs PC2") +
  stat_ellipse(level = 0.68) + 
  # geom_text(aes(label = label),
  #           data = adonis_label,
  #           vjust = "top", hjust = "right", color = "black") + 
  theme(axis.text = element_text(size = 11,color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5))
ggsave("PCA.png", path = Beta_diversity.path)

## PCoA
res <- ape::pcoa(distance_matrix, correction = "none", rn = NULL)
result <- res$values[, "Relative_eig"]
pro1 <- as.numeric(sprintf("%.3f",result[1])) * 100
pro2 <- as.numeric(sprintf("%.3f",result[2])) * 100
x <- res$vectors
sample_names <- rownames(x)
pc <- as.data.frame(res$vectors)
pc$names <- sample_names
legend_title <- ""
pc$group <- otu.group[rownames(pc),]

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
ggsave("PCoA.png", path = Beta_diversity.path)

## NMDS
nmds <- metaMDS(otu.rare, trace = F)
site.scores <- scores(nmds) %>% as.data.frame()
site.scores$group <- otu.group[rownames(site.scores),]

p <- ggplot(data = site.scores,aes(x=NMDS1,y=NMDS2)) + 
  geom_point(aes(shape=group,color=group,fill=group),size=3) + 
  labs(title = "NMDS - NMDS1 vs NMDS2") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept=0,linetype="dashed",color="blue")+geom_hline(yintercept = 0,linetype="dashed",color="blue") + 
  geom_text(x=0.2,y=-0.2,label=paste("stress = ",round(nmds$stress,3),sep = ""),color="blue",size=5)
ggsave("NMDS.png", path = Beta_diversity.path)

## Analysis of community structure difference between groups, anosim, MRPP, adonis
## anosim plot
otu.anosim <- anosim(distance_matrix, otu.group[rownames(distance_matrix),], permutations=999)
plt.pdf <- paste(Beta_diversity.path, "anosim.pdf", sep = "")
# pdf(plt.pdf)
# plot(otu.anosim, col=rainbow(group.count+1), main="Anosim", xlab="", ylab="")
# dev.off()

## Run together anosim, MRPP, adonis
## Just use NULL, ERROR: s--- rbind will as factor(matrix %>% as.dataframe)
# anosim.table <- adonis.table <- mrpp.table <- matrix()[-1,] %>% as.data.frame()
anosim.table <- adonis.table <- mrpp.table <-NULL
for (i in 1:(group.count-1)) {
  for (u in (i+1):group.count) {
    col1.temp <- paste(group.unique[i], "/", group.unique[u], sep=" ")
    otu.temp <- filter(otu.data, group == group.unique[i] | group == group.unique[u])
    
    ## anosim
    anosim.temp <- anosim(otu.temp[,-ncol(otu.temp)], otu.temp$group, permutations = 999, distance = 'bray')
    anosim.vector <- c(col1.temp,round(anosim.temp$statistic, digits = 5), anosim.temp$signif)
    anosim.table <- rbind(anosim.table, anosim.vector)
    
    ## adonis 
    adonis.temp <- adonis(otu.temp[,-ncol(otu.temp)]~otu.temp$group, permutations = 999, method = 'bray')
    adonis.temp.stat <- round(unlist(adonis.temp$aov.tab[1,], use.names = F), digits = 5) %>% as.character()
    Residuals <- round(unlist(adonis.temp$aov.tab[2,], use.names = F), digits = 5) %>% as.character()
    adonis.temp.out <- c()
    for (s in 1:6) {
      if (is.na(Residuals[s])) {
        adonis.temp.out <- c(adonis.temp.out,adonis.temp.stat[s])
      } else {
        adonis.temp.out <- c(adonis.temp.out, str_c(adonis.temp.stat[s],"(",Residuals[s],")"))
      }
    }
    adonis.temp.out <- c(col1.temp, adonis.temp.out)
    adonis.table <- rbind(adonis.table, adonis.temp.out)
    
    ## MRPP
    mrpp.temp <- mrpp(otu.temp[,-ncol(otu.temp)], otu.temp$group, permutations = 999, distance = 'bray')
    mrpp.temp.out <- round(c(mrpp.temp$A, mrpp.temp$delta, mrpp.temp$E.delta, mrpp.temp$Pvalue), digits = 5)
    mrpp.temp.out <- c(col1.temp, mrpp.temp.out)
    mrpp.table <- rbind(mrpp.table, mrpp.temp.out)
  }
}
colnames(anosim.table) <- c("Group", "ANOSIM-statistic-R", "Significance")
colnames(adonis.table) <- c("Group", dimnames(adonis.temp$aov.tab)[[2]])
colnames(mrpp.table) <- c("Group", "A", "observed-delta", "expected_delta", "Significance")

anosim.file <- str_c(Beta_diversity.path, "anosim.xls")
adonis.file <- str_c(Beta_diversity.path, "adonis.xls")
mrpp.file <- str_c(Beta_diversity.path, "mrpp.xls")
write.table(anosim.table, anosim.file, quote = F, sep="\t",row.names = F, col.names = T, fileEncoding = "UTF-8")
write.table(adonis.table, adonis.file, quote = F, sep="\t",row.names = F, col.names = T, fileEncoding = "UTF-8")
write.table(mrpp.table, mrpp.file, quote = F, sep="\t",row.names = F, col.names = T, fileEncoding = "UTF-8")

## different species between groups 
## ranksum & T-test
otu.id <- colnames(otu.rare)
otu.use <- as.data.frame(otu.rare)
otu.use$group <- otu.group[rownames(otu.use),]

if ( group.count < 3 ) {
  ## wilcoxon test
  ranksum.diff <- apply(otu.use[,-ncol(otu.use)], 2, my.wilcoxon)
  ## T test
  otu.differences <- apply(otu.use[,-ncol(otu.use)], 2, my.ttest)
} else {
  ## K-W test
  ranksum.diff <- apply(otu.use[,-ncol(otu.use)], 2, my.kruskal)
  ## aov differeces analysis
  otu.differences <- apply(otu.use[,-ncol(otu.use)], 2, my.aov)
}

func.sign <- function(ranksum.diff) {
  #~ [library] P adjust : stats or fdrtools   
  # ranksum.diff <- rbind(ranksum.diff, P.adjust=fdrtool(ranksum.diff[2,], statistic = 'pvalue', plot = F)$qval)
  ranksum.diff <- rbind(ranksum.diff, P.adjust=p.adjust(ranksum.diff[2,], method = 'fdr'))
  ## p.adjust("fdr") < 0.1, significant differences
  # otu.significant <- subset(as.data.frame(t(ranksum.diff)), P.adjust < 0.1)
  otu.significant <- as.data.frame(t(ranksum.diff))
  ## taxonomy id --> scientific name
  otu.significant <- cbind(Taxa=rownames(otu.significant), otu.significant)
  otu.significant <- otu.significant[order(otu.significant[,3]),]
  colnames(otu.significant)[3] <- "P.value"
  return(otu.significant)
}

otu.significant <- func.sign(ranksum.diff)
otu.significant2 <- func.sign(otu.differences)

otu.significant.file <- str_c(Beta_diversity.path, "ranksum_taxa_significant.xls")
write.table(otu.significant, otu.significant.file, quote = F, sep="\t",row.names = F, fileEncoding = "UTF-8")
otu.significant2.file <- str_c(Beta_diversity.path, "differences_taxa_significant.xls")
write.table(otu.significant2, otu.significant2.file, quote = F, sep="\t",row.names = F, fileEncoding = "UTF-8")

## metagenomeseq
#@ ERROR loadMeta: sample names can't be digital 
# otu.meta <- loadMeta(str_c(Intermediate.path, "meta.taxa_counts.csv"))
otu.meta_counts <- read.delim(str_c(Intermediate.path, "meta.taxa_counts.csv"), row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% as.data.frame()
otu.meta_taxa <- data.frame(taxa=rownames(otu.meta_counts))
otu.meta <- list(counts=otu.meta_counts, taxa=otu.meta_taxa)
taxa <- read.delim(str_c(Intermediate.path, "meta.taxa_taxanomy.csv"), stringsAsFactors = F)
#@ group.txt header can't startswith '#'!!!
clin <- loadPhenoData(meta.group_file, tran = T)
clin$siteSampled <- "intestine"  ## Placeholder, prevent data.frame --> character
clin$group <- gsub('group', 'GROUP', clin$group)
clin.unique <- unique(clin$group)
rareFeatures_threshold <- length(clin$group) / 66 * 10

ord <- match(colnames(otu.meta$counts), rownames(clin))
clin <- clin[ord, ]
phenotypeData <- AnnotatedDataFrame(clin)
OTUdata <- AnnotatedDataFrame(taxa)
obj <- newMRexperiment(otu.meta$counts, phenoData = phenotypeData, featureData = OTUdata)
rareFeatures <- which(rowSums(MRcounts(obj) > 0) < rareFeatures_threshold)
if ( length(rareFeatures) != 0 ) {obj <- obj[-rareFeatures, ]}  ## df[-0,]  ----  remove all rows (all data)
obj <- filterData(obj, present = ncol(otu.meta$counts)*0.375, depth = 1)
objp <- cumNormStat(obj, pFlag = F)
obj <- cumNorm(obj, p=objp)

if (group.count == 2) {
  ## 2 groups 
  pd <- pData(obj)
  mod <- model.matrix(~group, data = pd)
  objresl <- fitFeatureModel(obj, mod)
  metaseq.top20 <- MRcoefs(objresl, number = 20)
  metaseq.top20_simple <- metaseq.top20[,1:4]
} else {
  ## multiple groups 
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

metaseq.result.file <- str_c(Beta_diversity.path, "metagenomeSeq.taxa_significant.xls")
metaseq.simple_result.file <- str_c(Beta_diversity.path, "metagenomeSeq.taxa_significant_simple.xls")
write.table(metaseq.top20, file = metaseq.result.file, sep = '\t', row.names = F, fileEncoding = 'UTF-8')
write.table(metaseq.top20_simple, file = metaseq.simple_result.file, sep = '\t', row.names = F, fileEncoding = 'UTF-8')
