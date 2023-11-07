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

sum(keep)/nrow(COUNTS) # 0.8129889 % N transcripts

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

RES.P <- RES %>% filter( padj < 0.05 & abs(log2FoldChange) > 2) 

RES.P %>% dplyr::count(sampleB)

UPSETDF <- RES.P %>% 
  # mutate(SIGN = sign(log2FoldChange)) %>%
  filter(log2FoldChange < 0 ) %>% # ONLY UP-EXPRESSED IN EXPERIMENTAL CNTRST (i.e. DOWN-EXP. IN Ref)
  group_by(Name) %>%
  summarise(across(sampleB, .fns = list), n = n()) 

library(ggupset)

UPSETDF %>%
  ggplot(aes(x = sampleB)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  # axis_combmatrix(levels = c("24 HPF", "110 HPF")) +
  labs(x = '', y = 'Number of transcripts (up-expressed)') +
  # scale_color_manual("", values = col) +
  # scale_fill_manual("", values =  col) +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p1

p1 <- p1 + theme(legend.position = "top",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

p1

UPSETDF %>% filter(n == 1) %>% dplyr::count(sampleB)
UPSETDF %>% filter(n == 1) %>% unnest(sampleB)
