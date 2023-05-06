#!/usr/bin/env Rscript

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
rgi.file <- paste(work.dir, "rgi_single.txt", sep = "/")
kraken.file <- paste(results.dir, "4.species_annotation", sample.num, 
                     "gene_annotation_by_kraken2.txt", sep = "/")
#~ check
stopifnot(dir.exists(work.dir),
          file.exists(rgi.file),
          file.exists(unigenes_id.file)
          )
print(paste("INFO", date(), paste0("function_annotation_card.r Work Path: ", work.dir), sep = " - "))
##############################################################################
## library
suppressMessages({
  library(tidyverse)
  library(ggsci)
  library(reshape)
  library(circlize)
})

# configure
my.palette <- pal_npg()(10)
my.pal27 <- c(pal_npg("nrc", alpha = 1)(9), pal_npg("nrc", alpha = 0.7)(9), pal_npg("nrc", alpha = 0.4)(9))

# function
#~ write table to workpath
my.write_table <- function(x, file_name, ...) {
  write.table(x = x, file = paste(work.dir, file_name, sep = "/"), quote = F, row.names = F, 
              sep = "\t", fileEncoding = "UTF-8", ...)
}
##############################################################################

## annotate
print(paste("INFO", date(), "Start annotate", sep = " - "))
#~ unigenes number
cmd.tmp <- paste0("wc -l ", unigenes_id.file)
unigenes.number <- strsplit(system(cmd.tmp, intern = T), " ")[[1]][1] %>% 
                    as.numeric()
## get unigenes phylum
kraken.table <- read.delim(kraken.file, header = T, sep = "\t", check.names = F, 
                           stringsAsFactors = F, row.names = 1)
phylum.vector <- kraken.table$phylum
names(phylum.vector) <- rownames(kraken.table)

## rgi result
#~ V1 ORF_ID; V9 ARO; V16 Resistance Mechanism
rgi.table <- read.delim(rgi.file, header = F, sep = "\t", check.names = F, stringsAsFactors = F)[,c(1,9,16)]
annotated_on_card.number <- dim(rgi.table)[1]
rgi.table$phylum <- phylum.vector[rgi.table$V1]
out.table <- rgi.table
colnames(out.table)[1:3] <- c("ORF_ID","ARO","Resistance_Mechanism")
my.write_table(out.table, "aro_rgi_simple.xls")

## ARO
aro.table <- rgi.table %>%
  select(V1,V9) %>%
  group_by(V9) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  mutate(ppm = round(count / unigenes.number * 10^6, digits = 2))
annotated_on_aro.number <- dim(aro.table)[1]
colnames(aro.table)[1] <- "ARO"
out.table <- aro.table[,-2]
my.write_table(out.table, "aro.xls")

## Resistance Mechanism
rm.table <- rgi.table %>% 
  select(V16, phylum) %>% 
  filter(startsWith(V16, "antibiotic") & phylum != "None")
#~ select first resistance mechanism; replace "antibiotic"
rm.table$V16 <- sapply(rm.table$V16, function(x) {
  str_replace(strsplit(x, ";")[[1]][1], pattern = "antibiotic ", replacement = "")
  })
rm.table <- rm.table %>% 
  group_by(V16, phylum) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))
## Resistance Mechanism chord diagram
#~ > 4 phylum, use Others replace
if (n_distinct(rm.table$phylum) > 4) {
  rm.table_cast <- cast(rm.table, phylum~V16, fill = 0)
  rm.table_cast <- rm.table_cast %>% arrange(desc(rowSums(rm.table_cast)))
  rownames(rm.table_cast) <- rm.table_cast$phylum
  rm.table_cast <- rm.table_cast[,-1]
  rm.table_cast <- rbind(rm.table_cast[1:4,], Others=colSums(rm.table_cast[5:nrow(rm.table_cast),]))
  rm.table_cast <- cbind(phylum=rownames(rm.table_cast), rm.table_cast)
  fig.table <- melt(rm.table_cast, id ="phylum", drop.margins = T)%>% 
    filter(value != 0)
} else {
  fig.table <- rm.table
}
color.number <- n_distinct(fig.table[,1]) + n_distinct(fig.table[,2])
#~ plot
png(paste(work.dir, "aro_phylum_resistance_mechanism_chord_diagram.png", sep = '/'),
units="in", width=9, height=9, res=500)
# par(cex = 1.2)
chordDiagram(fig.table, scale = T, annotationTrack = c("grid", "name"), 
             grid.col = pal_npg()(color.number), col = pal_npg()(color.number))
circos.clear()
dev.off()

## card function annotation stat
col1.tmp <- c("Gene catalogue:", "Annotated on CARD:", "Annotated AROs:")
col2.tmp <- c(unigenes.number, annotated_on_card.number, annotated_on_aro.number)
out.table <- data.frame(col1.tmp, col2.tmp)
my.write_table(out.table, "function_annotation_stats_card.xls", col.names = F)
