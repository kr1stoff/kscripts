#!/usr/bin/env Rscript


## get work path
## get work path
print_Usage <- function() {
  cat("PROG <results_dir> <sample.number>\n")
  cat("tips:\n\tresults_dir ==> ./results\n\tsample.number ==> SRR11575973\n")
  quit("no", status = 1)
}
args <- commandArgs(T)
if (!length(args) == 2) {print_Usage()}
results.dir <- args[1]           # /nfs-test/mxf/Project/2.metagenomics/202010_metapipe/results
sample.num <- args[2]            # SRR11575973
work.dir <- paste(results.dir, "5.function_annotation", sample.num, sep = "/")
unigenes_id.file <- paste(results.dir, "3.gene_prediction/bowtie2_out", sample.num, 
                          "unigenes_id.protein.txt", sep = "/")
#~ check
stopifnot(length(args) == 2,
          dir.exists(work.dir),
          file.exists(paste(work.dir, 'out.emapper_simgle.annotations', sep = '/')),
          file.exists(unigenes_id.file)
          )

print(paste("INFO", date(), paste0("function_annotation.r Work Path: ", work.dir), sep = ' - '))
# db path
kegg_db_path <- '/nfs-test/mxf/Database/KEGG/kegg_BRITE_table.txt'
nog_db_path <- '/nfs-test/mxf/Database/NCBI/COG/Functional_Category.txt'
cazy_db_path <- '/nfs-test/mxf/Database/CAZy/cazy.ref.txt' 
go_db_path <- '/nfs-test/mxf/Database/GO/go_comparison.txt'

##############################################################################
## library
suppressMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggsci)
  library(scales)
})

# configure
my.palette <- pal_npg()(10)
my.pal27 <- c(pal_npg("nrc", alpha = 1)(9), pal_npg("nrc", alpha = 0.7)(9), pal_npg("nrc", alpha = 0.4)(9))

# function
#~ write table to workpath
my.write_table <- function(x, file_name, ...) {
  write.table(x = x, file = paste(work.dir, file_name, sep = '/'), quote = F, row.names = F, sep = '\t', fileEncoding = 'UTF-8', ...)
}
#~ emapper indices(go,ko,cazy..) extract, calculate number of matched genes
my.emapper.stat <- function(colnum, index, rp_mark=NULL) {
  intermedia <- emapper.result[,c(1,colnum)]
  colnames(intermedia) <- c("V1", "V2")
  intermedia <- intermedia %>%
    group_by(V1) %>%
    filter(V2 != "")
  if (is.null(rp_mark)) {
    intermedia <- transmute(intermedia, func = str_split(V2, ',', simplify = T)[1,1])
  } else {
    intermedia <- transmute(intermedia, func = str_replace(str_split(V2, ',', simplify = T)[1,1], rp_mark, ""))
  }
  emmaper.func <- intermedia %>%
    group_by(func) %>%
    summarise(number_of_matched_genes = n()) %>%
    arrange(desc(number_of_matched_genes))
  colnames(emmaper.func)[1] <- index
  return(emmaper.func)
}

##############################################################################

## emapper results
emapper.file <- paste(work.dir, 'out.emapper_simgle.annotations', sep = '/')
emapper.result <- read.delim(emapper.file, header = F, sep = '\t', check.names = F, stringsAsFactors = F)
#~ unigenes number
# number.all_genes <- dim(emapper.result)[1]
cmd.tmp <- paste0("wc -l ", unigenes_id.file)
number.all_genes <- strsplit(system(cmd.tmp, intern = T), " ")[[1]][1] %>% 
  as.numeric()

## GO annotation
emapper.go <- my.emapper.stat(7, index = 'go')
number.go_genes <- sum(emapper.go$number_of_matched_genes)
my.write_table(emapper.go, "go.xls")

go.terms <- emapper.go$go
gene <- mapIds(org.Hs.eg.db, keys = go.terms, column = 'ENTREZID', 'SYMBOL', keytype = "GO")
gene <- gene[which(!is.na(gene))]

run_gf <- function(gf) {
  GF.params <- enrichGO(gene = gene,    
                        OrgDb = org.Hs.eg.db,    
                        ont = gf,    
                        pAdjustMethod = "fdr",    
                        pvalueCutoff = 0.01,    
                        qvalueCutoff  = 0.05)
  df_temp <- head(GF.params@result[,c("Description", "GeneRatio", "p.adjust", "Count")],10) %>% mutate(ONTOLOGY=gf)
  return(df_temp)
}
gf.list <- lapply(c("BP","CC","MF"), run_gf)
#~ goALL dataframe for plot
goAll <- data.frame()
for (gfitem in gf.list) {goAll <- rbind(goAll, gfitem)}
goAll$ONTOLOGY_Description <- paste0(goAll$Description, "[", goAll$ONTOLOGY, ']')
#~ string to decimal
calc_fraction <- function(x) {eval(parse(text = x))}
goAll$GeneRatio <- sapply(goAll$GeneRatio, calc_fraction)
my.write_table(goAll, "go.enrich.xls")
#~ GF barplot
p <- ggplot(goAll) +
  geom_bar(aes(x = Description, y = -log10(p.adjust), fill = ONTOLOGY), stat = 'identity') +
  coord_flip() +
  labs(x=NULL, title = "GO Function Catalogue annotation") + 
  scale_x_discrete(limits=goAll$Description) +
  scale_fill_manual(values = my.palette) +
  theme_bw()
ggsave(paste(work.dir, "go.barplot.png", sep = '/'), p, width = 12, height = 9)
#~ GF point
p <- ggplot(goAll, aes(x = GeneRatio, y = ONTOLOGY_Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_y_discrete(limits = goAll$ONTOLOGY_Description) +
  # ggtitle("GO enrichment") +
  scale_color_gradient(low = my.palette[1], high = my.palette[3]) +
  labs(y = NULL, title = "GO Enrichment Analysis") +
  theme_bw()
ggsave(paste(work.dir, "go.pointplot.png", sep = '/'), p, width = 12, height = 9)

## KEGG annotation
#~ ko
kegg.ko <- my.emapper.stat(9, index = 'ko', rp_mark = "ko:")
number.ko_genes <- sum(kegg.ko$number_of_matched_genes)
my.write_table(kegg.ko, "kegg.ko.xls")
#~ ec
kegg.ec <- my.emapper.stat(8, index = 'ec')
number.ec_gens <- sum(kegg.ec$number_of_matched_genes)
my.write_table(kegg.ec, "kegg.ec.xls")
#~ module
kegg.module <- my.emapper.stat(11, index = 'module')
number.module_genes <- sum(kegg.module$number_of_matched_genes)
my.write_table(kegg.module, "kegg.module.xls")
#~ kegg ko00001 BRITE table
kegg.brite <- read.delim(kegg_db_path, sep = '\t', header = T, check.names = F, stringsAsFactors = F, fileEncoding = "UTF-8", colClasses = "character", quote = "")[,1:6] %>%
  unique()
#~~ pathway summary
kegg.pathway <- my.emapper.stat(10, index = 'LevelC_BRITE_ID', rp_mark = 'ko')
number.pathway_genes <- sum(kegg.pathway$number_of_matched_genes)
kegg.summary <- merge(kegg.brite, kegg.pathway, by = 'LevelC_BRITE_ID')
kegg.summary <- dplyr::select(kegg.summary, colnames(kegg.brite), everything())
kegg.summary_6lv1 <- filter(kegg.summary, !LevelA_BRITE_NAME %in% c("Brite Hierarchies", "Not Included in Pathway or Brite"))
#~ pathway lv3
kegg.lv3 <- dplyr::select(kegg.summary_6lv1, LevelC_BRITE_NAME, number_of_matched_genes) %>%
  group_by(LevelC_BRITE_NAME) %>%
  arrange(desc(number_of_matched_genes))
my.write_table(kegg.lv3, "kegg.pathway_lv3.xls")
#~ pathway lv2
kegg.lv2 <- dplyr::select(kegg.summary_6lv1, LevelB_BRITE_NAME, number_of_matched_genes) %>%
  group_by(LevelB_BRITE_NAME) %>%
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>%
  arrange(desc(number_of_matched_genes))
my.write_table(kegg.lv2, "kegg.pathway_lv2.xls")
#~ pathway lv1
kegg.lv1 <- dplyr::select(kegg.summary_6lv1, LevelA_BRITE_NAME, number_of_matched_genes) %>%
  group_by(LevelA_BRITE_NAME) %>%
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>%
  mutate(fraction_of_matched_genes = number_of_matched_genes / number.all_genes) %>%
  arrange(desc(number_of_matched_genes))
my.write_table(kegg.lv1, "kegg.pathway_lv1.xls")
#~ group barplot
df.fig <- dplyr::select(kegg.summary_6lv1, LevelA_BRITE_NAME, LevelB_BRITE_NAME, number_of_matched_genes) %>%
  group_by(LevelA_BRITE_NAME, LevelB_BRITE_NAME) %>%
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>% 
  arrange(LevelA_BRITE_NAME)
df.fig$LevelB_BRITE_NAME <- factor(df.fig$LevelB_BRITE_NAME, levels = rev(df.fig$LevelB_BRITE_NAME))

p <- ggplot(df.fig, aes(x = number_of_matched_genes, y = LevelB_BRITE_NAME, fill = LevelA_BRITE_NAME)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  labs(title = "KEGG pathway annotation",
       x = "Number of matched genes",
       y = NULL,
       fill = NULL) + 
  scale_fill_manual(values = my.pal27)
  # geom_text(data = df.fig, aes(label = number_of_matched_genes, hjust = 'right'), size = 3) # mute text at bar
ggsave(filename = str_c(work.dir,"kegg_pathway.png", sep = '/'), plot = p, dpi = 300, width = 12, height = 9)

## NOG annotation
#~ OGs
eggnog.og <- my.emapper.stat(19, index = 'og', rp_mark = '@.*')
number.og_genes <- sum(eggnog.og$number_of_matched_genes)
temp.add_nog_header <- function(x) {
  if (grepl('COG', x)){
    return(x)
  } else {
    return(paste0('ENOG50', x))
  }
}
eggnog.og$og <- sapply(eggnog.og$og, temp.add_nog_header)
my.write_table(eggnog.og, "eggnog.og.xls")
#~ NOG funtion
eggnog.func <- my.emapper.stat(21, index = 'nog_func')
number.func_genes <- sum(eggnog.func$number_of_matched_genes)
eggnog.func$nog_func <- substr(eggnog.func$nog_func, 1, 1) 
eggnog.func <- group_by(eggnog.func, nog_func) %>%
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>%
  arrange(desc(number_of_matched_genes))
nog.db <-  read.table(nog_db_path, sep = '\t', header = F, col.names = c("abbr_lv2","func_lv2_name","func_lv1_name"), check.names = F, stringsAsFactors = F, fileEncoding = "UTF-8")
nog.summary <- merge(nog.db, eggnog.func, by.x = 'abbr_lv2', by.y = 'nog_func') %>%
  mutate(fraction_of_matched_genes = number_of_matched_genes / number.all_genes) %>%
  mutate(plot_legend = paste(abbr_lv2,func_lv2_name,sep = ' : ')) %>%
  arrange(desc(number_of_matched_genes))
# my.write_table(nog.summary[,-ncol(nog.summary)], "eggnog.function.xls")
#~ eggnog lv2
nog.lv2 <-  nog.summary[,c('plot_legend', 'number_of_matched_genes', 'fraction_of_matched_genes')]
colnames(nog.lv2)[1] <- 'func_lv2_name'
my.write_table(nog.lv2, "eggnog.function_lv2.xls")
#~ eggnog lv1
nog.lv1 <- group_by(nog.summary, func_lv1_name) %>%
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>%
  arrange(desc(number_of_matched_genes))
my.write_table(nog.lv1, "eggnog.function_lv1.xls")
#~ plot
p <- ggplot(nog.summary, mapping = aes(x = abbr_lv2, y = number_of_matched_genes, fill = plot_legend)) + 
  geom_bar(stat = 'identity') +
  theme_bw() + 
  labs(title = "eggNOG Function Category annotation",
       x = NULL,
       y = "Count of matched genes",
       fill = NULL) + 
  scale_fill_manual(values = my.pal27)
  # geom_text(mapping = aes(label = number_of_matched_genes, vjust = "bottom"), size = 3)
ggsave(str_c(work.dir, "eggnog_func.png", sep = '/'), plot = p, height = 9, width = 14)

## CAZy
cazy.ref <- read.table(cazy_db_path, sep = '\t', stringsAsFactors = F, check.names = F, header = T)
cazy <- my.emapper.stat(16, index = 'cazy')
number.cazy_genes <- sum(cazy$number_of_matched_genes)
my.write_table(cazy, "cazy_func_lv2.xls")
cazy.summary <- merge(cazy.ref, cazy, by.x = 'level_2', by.y = 'cazy')
#~ level 1
level1.table <- group_by(cazy.summary, level_1_abbr, level_1_name) %>% 
  summarise(number_of_matched_genes = sum(number_of_matched_genes)) %>%
  mutate(fraction_of_matched_genes = number_of_matched_genes / number.all_genes) %>%
  arrange(desc(number_of_matched_genes))
cazy.plot <- level1.table %>% mutate(plot_legend = paste(level_1_abbr, level_1_name, sep = ' : '))
level1.table_out <- cazy.plot[,c('plot_legend','number_of_matched_genes','fraction_of_matched_genes')]
colnames(level1.table_out)[1] <- 'cazy_lv1'
my.write_table(level1.table_out, "cazy_func_lv1.xls")
#~ CAZy level 1 plot
p <- ggplot(cazy.plot, mapping = aes(x = level_1_abbr, y = number_of_matched_genes, fill = plot_legend)) +
  geom_bar(stat = 'identity') +
  theme_bw() + 
  labs(title = "CAZy Function annotation",
       x = NULL,
       y = "Count of matched genes",
       fill = NULL) + 
  scale_fill_manual(values = my.pal27)
  # geom_text(mapping = aes(label = number_of_matched_genes, vjust = "bottom"), size = 3)
ggsave(str_c(work.dir, "cazy_func.png", sep = '/'), plot = p, height = 7, width = 7.5)

## annotation stat
my.stat_format <- function(number_of_genes) {
  return(paste0(number_of_genes, "(", percent(number_of_genes/number.all_genes, accuracy = 0.01), ")"))
}
number.indices <- c(number.go_genes, number.ko_genes, number.ec_gens, number.module_genes, number.pathway_genes, number.og_genes, number.func_genes, number.cazy_genes)
stat.col2 <- sapply(number.indices, my.stat_format)
stat.col2 <- c(number.all_genes, stat.col2)
stat.col1 <- c("Gene catalogue:", "Annotated on GO:", "Annotated on KO:", "Annotated on EC:", 
               "Annotated on Module:", "Annotated on pathway:", "Annotated on OG:", 
               "Annotated on eggNOG:", "Annotated on CAZy:")
stat.table <- data.frame(stat.col1, stat.col2)
my.write_table(stat.table, 'function_annotation_stats.xls', col.names = F)
