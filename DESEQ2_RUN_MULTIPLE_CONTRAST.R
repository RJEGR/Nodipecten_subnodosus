#
# RUN MULTIPLE CONTRAST D.E ANALYSIS
# RICARDO GOMEZ-REYES
# 2023

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MANO_DELEON/"

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

count_f <- list.files(path = path, pattern = "count.isoform.counts.matrix", full.names = T)

MTD_f <- list.files(path = path, pattern = "metadata_Pau.tsv$", full.names = T)


.colData <- read_tsv(MTD_f)

COUNTS <- read_tsv(count_f)

colnames(COUNTS)[1] <- "ID"

rowNames <- COUNTS$ID

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

# .COUNTS <- COUNTS

# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep) # N transcripts

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)

# 2) Run Multiple Contrast Comparison =====

library(DESeq2)

CONTRAST <- .colData %>% dplyr::select(starts_with("CONTRAST")) %>% names()

run_contrast_DE <- function(COUNTS, colData, CONTRAST = NULL, ref = NULL) {
  
  # CONTRAST: Column in colData with character vector of Design
  
  names(colData)[1] <- "LIBRARY_ID"
  
  names(colData)[names(colData) %in% CONTRAST] <- "Design"
  
  colData <- colData %>% drop_na(Design)
  
  # any(colnames(COUNTS) == colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  keep <- colnames(COUNTS) %in% colData$LIBRARY_ID 
  
  COUNTS <- COUNTS[,keep]
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = ref))
  
  require(DESeq2)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = COUNTS,
    colData = colData,
    design = ~ Design )
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  dds <- estimateDispersions(dds)
  
  dds <- nbinomWaldTest(dds)
  
  # return(dds)
  
  contrast <- levels(colData(dds)$Design)
  
  res <- get_res(dds, contrast)
  
  return(res)
}

# run_contrast_DE(COUNTS, .colData, CONTRAST = CONTRAST[1], ref = "Ref")

any(colnames(COUNTS) == .colData$Library_ID)

out <- list()

for (j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  
  out [[i]] <- run_contrast_DE(COUNTS, .colData, CONTRAST = CONTRAST[i], ref = "Ref") %>% 
    dplyr::mutate(CONTRAST = CONTRAST[i])
}

do.call(rbind, out) -> RES

# AFTER THEN, EXPORT RESULTS FILE AND CONTINUE W/ GENE ONTOLOGY ENRICHMENT ANALYSIS
