# Evaluate multiple contrast of the experimetnal design, the subset count-matrix ALL degs and calculate lm to estimate effect size of contrasted groups. This is an step to assess what contrast are pivotal to describe

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

dir <- "~/Documents/MANO_DELEON/DATOS_Paulina_dir/"

subdir <- ""

f <- "isoforms.counts_experimental.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)

library(tidyverse)

dim(datExpr <- read.delim(f, row.names = 1))

datExpr <- round(datExpr)

.colData <- list.files(dir, pattern = "metadata_experimental.tsv", full.names = T) 

.colData <- read_tsv(.colData) %>% filter(Condition != "Reg")

table(.colData$Site, .colData$Condition)

.colData <- 
  .colData %>% mutate(sample_pairs = paste(Site, Condition, sep = "_")) %>%
  mutate_if(is.character, as.factor)


run_DESEQ <- function(count, contrast_col, colData) {
  
  names(colData)[names(colData) %in% "design"] <- "drop_duplicate_col"
  
  colData <- colData %>% dplyr::rename("design" = contrast_col)
  
  g <- levels(colData$design)
  
  create_pairs <- function(vec) {
    combn(vec, 2, simplify = FALSE)
  }
  
  
  combinations_vector <- unlist(lapply( create_pairs(g), paste, collapse = "-"))
  
  # Print the combinations vector
  print(combinations_vector)
  
  gstr <- strsplit(combinations_vector, "-")
  
  
  run_analysis <- function(count, gstr, colData) {
    
    cat("\nContrasting: ", gstr, "\n")
    
    colData <- colData %>% filter(design %in% gstr)
    
    sam_gstr <- structure(
      as.character(colData$design), 
      names = as.character(colData$Library_ID))
    
    sam_gstr <- sam_gstr[sam_gstr %in% gstr]
    
    fl <- as.factor(sam_gstr)
    
    colData <- colData %>% mutate(design = relevel(design, ref = levels(fl)[1]))
    
    
    # design <- model.matrix(~fl)
    
    # filter samples from the contrast
    
    keep_cols <- colnames(count) %in% names(sam_gstr)
    
    count <- round(count[,keep_cols])
    
    by_count <- 1; by_freq <- 1
    
    keep_genes <- rowSums(count > by_count) >= by_freq
    
    cat("\nUsing N genes: ", sum(keep_genes),"\n")
    
    cat("\nUsing samples: ", colnames(count),"\n")
    
    # dim(count <- count[keep_genes, keep_cols ])
    
    colnames(count) <- sam_gstr
    
    count <- as(count, 'matrix')
    
    count <- round(count)
    require(DESeq2)
    
    ddsFullCountTable <- DESeqDataSetFromMatrix(
      countData = count,
      colData = colData,
      design = ~ design)
    
    dds <- DESeq(ddsFullCountTable)
    
    result_table = results(dds)
    
    
    sampleA <- levels(fl)[1]
    sampleB <- levels(fl)[2]
    
    result_table <- data.frame(sampleA,
      sampleB,
      result_table)
    
    
    reorder_cols <- c( "ids" ,"sampleA", "sampleB", "log2FoldChange",  "pvalue", "padj")
    
    result_table %>%
      as_tibble(rownames = 'ids') %>%
      mutate_at(vars(!matches("ids|sample|PValue|FDR")),
        round ,digits = 2) %>%
      select_at(vars(all_of(reorder_cols)))
  }
  
  OUT <- lapply(gstr, function(x) run_analysis(count, gstr = x, colData))
  
  do.call(rbind, OUT) %>% filter(padj < 0.05)
  
}


OUT <- run_DESEQ(datExpr, contrast_col = "sample_pairs", colData = .colData)

write_rds(OUT, file = paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds"))

# Multiple contrast by Site (Evals si localidad tiene un efecto en el estres)


OUT_condition <- run_DESEQ(datExpr, contrast_col = "Condition", colData = .colData)

OUT_condition %>% dplyr::count(sampleA, sampleB) 

write_rds(OUT_condition, file = paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_condition.rds"))

# .colData %>% dplyr::count(Site)

# OUT_site <- run_DESEQ(datExpr, contrast_col = "Site", colData = .colData)

# write_rds(OUT_site, file = paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_site.rds"))

