

# EDA of 

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

Manifest <- .colData %>% mutate(design = paste(Site, Condition, sep = "_")) %>%
  mutate_if(is.character, as.factor)


DEGS <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_condition.rds")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

DEGS <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds")) %>%
  filter(sampleA != "BLA_Reg" & !sampleB  %in% c("BLA_Reg", "LOL_Reg")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)


read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_condition.rds")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05) %>%
  write_tsv(file = file.path(dir, subdir, "/p05_DESEQ_condition_analyis.tsv"))

read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds")) %>%
  filter(sampleA != "BLA_Reg" & !sampleB  %in% c("BLA_Reg", "LOL_Reg")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05) %>%
  write_tsv(file = file.path(dir, subdir, "/p05_DESEQ_multiple_contrast_analysis.tsv"))


# rbind(DEGS_condition, DEGS) 

# View by condition


DataVizdf <- DEGS %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2 ) %>%
  mutate(facet = ifelse( sign(log2FoldChange) == 1, "up in sampleA", "up in sampleB")) %>%
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  dplyr::count(sampleX, sort = T) 

# Estimar La frecuencia de degs, en cada contraste, resumindo a los factores del disenio experimental
degs_df <- DataVizdf %>% 
  group_by(sampleX) %>%
  summarise(degs = sum(n), n_contrast = n()) %>% arrange(desc(degs))


UPSETDF <- DEGS %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2 ) %>%
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  # separate(x, into = c("Site", "Condition"), sep = "_", remove = F) %>%
  # count(x, Site, Condition)
  # group_by(ids, Site, Condition) %>%
  distinct(ids, sampleX) %>%
  group_by(ids) %>%
  summarise(across(sampleX, .fns = list), n = n()) 


UPSETDF %>% filter(n == 1) %>% unnest(sampleX) %>% 
  dplyr::count(sampleX, sort = T) %>%
  dplyr::rename("unique_degs" = "n") %>%
  left_join(degs_df) %>% view()


library(ggVennDiagram)

DF <- UPSETDF %>%  unnest(sampleX)

DF %>% dplyr::count(sampleX)

recode_design <- c(
  # = "Differences in populations (Basal site)"
  "Bas" = "Basal",
  "Cao" = "Global acclimatation\nto Cao challenge",
  "Cte" = "Global acclimatation\nto Cte challenge",
  "BLA_Bas" = "Basal",
  "BLA_Cao" = "Specific acclim.\nto chaotic challenge (BLA)",
  "BLA_Cte" = "Specific acclim.\nto cte challenge (BLA)",
  "LOL_Bas" = "Basal",
  "LOL_Cao" = "Specific acclim.\nto chaotic challenge (LOL)",
  "LOL_Cte" = "Specific acclim.\nto cte challenge (LOL)")


# DF <- DF %>% 
#   dplyr::mutate(sampleX = dplyr::recode(sampleX, !!!recode_design)) %>% 
#   filter(sampleX != "Basal")

gene2ven <- split(DF$ids, DF$sampleX)

# gene2ven <- split(strsplit(UPSETDF$ids, "") , DF$Design)

gene2ven <- lapply(gene2ven, unlist)

str(gene2ven)

ggVennDiagram(gene2ven,label_font = "GillSans", label_size = 7,
  relative_height = 1,relative_width = 2) + 
  scale_fill_gradient(low="grey90",high = "red",)

ggVennDiagram(gene2ven, force_upset = T, order.intersect.by = "size")

