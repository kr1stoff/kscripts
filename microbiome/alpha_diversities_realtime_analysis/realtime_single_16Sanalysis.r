#!/usr/bin/env R

library(getopt)
spec <- matrix(
  c("input", "i", 1, "character", "Input file, like OTU/Species/Genus abundance table.",
    "outpath", "o", 2, "character", "Result files path. [default: ./16S_result_realtime]",
    "help", "h", 0, "logical", "Display help infomation."
  ),
  byrow = T, ncol = 5)
opt <- getopt(spec)
if ( !is.null(opt$help) || is.null(opt$input)) {
  cat(getopt(spec, usage = T))
  q(status = 1)
}
if ( is.null(opt$outpath) ) { opt$outpath <- "./16S_result_realtime" }
cat("\nInput ", opt$input, 
    "\nOutpath ", opt$outpath, "\nConfiguration is successful!\n")

suppressMessages({
library(vegan)
library(tidyverse)
})
############################# MY FUNCTION ######################################
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)	# richness index 
  else if (method == 'chao1') result <- estimateR(x)[2, ]	#Chao1 index
  else if (method == 'ace') result <- estimateR(x)[4, ]	#ACE index
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)	#Shannon index
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')	#Gini-Simpson index
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou uniformity
  else if (method == 'goods_coverage') result <- 1 - rowSums(x == 1) / rowSums(x)	#good's coverage
  else if (method == 'pd' & !is.null(tree)) {	#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  names(result) <- NULL
  result
}
func.apply1 <- function(method) {alpha_index(x=otu, method = method)}
create.diretory <- function(path){
  dir.create(str_c(path, 'result_old', sep = '/'), recursive = T)
  dir.create(str_c(path, 'result', sep = '/'), recursive = T)
}
func.create_result <- function(otu.alpha) {
  dir.create(fold.new, recursive = T)
  write.table(otu.alpha, file = outputfile.new, sep = '\t', row.names = F, quote = F, fileEncoding = 'UTF-8')
  cycles <- 1:nrow(otu.alpha)
  ## plot alpha index
  for (i in 1:length(alpha.index_list)) {
    alpha.index <- alpha.index_list[i]
    fig.name <- str_c(opt$outpath, '/result/', alpha.index, ".png")
    png(fig.name)
    plot(cycles, otu.alpha[,i], type = 'o', pch = 20, ylab = '', main = alpha.index)
    dev.off()
  }
}
################################################################################

otu <- read.delim(opt$input, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% t()
alpha.index_list <- c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage")
otu.alpha <- sapply(alpha.index_list, func.apply1) %>% t() %>% as.data.frame()

## collate output files
# if ( !dir.exists(opt$outpath) ) {create.diretory(opt$outpath)} 
fold.old <- str_c(opt$outpath, 'result_old', sep = '/')
fold.new <- str_c(opt$outpath, 'result', sep = '/')
outputfile.old <- str_c(fold.old, 'taxa_alpha.txt', sep = '/')
outputfile.new <- str_c(fold.new, 'taxa_alpha.txt', sep = '/')
## old directery 
if (!dir.exists(fold.new)) {
  func.create_result(otu.alpha)
} else {
  alpha.old <- read.table(outputfile.new, check.names = F, stringsAsFactors = F, header = T)
  otu.alpha <- rbind(alpha.old, otu.alpha)
  if (dir.exists(fold.old)) {system(str_c("rm -r ", fold.old))}
  system(str_c("mv", fold.new, fold.old, sep = ' '))
  func.create_result(otu.alpha)
}
