#
# RUN GO ENRICHMENT ANALYSIS 
# USING TOPGO AND REVIGO 
# RICARDO GOMEZ-REYES
# 2023

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MANO_DELEON/"

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

RES.P <- read_rds(paste0(path, "/DESEQ2RES.rds")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 2) 

QUERIES <- RES.P %>% 
  # mutate(SIGN = sign(log2FoldChange)) %>%
  filter(log2FoldChange < 0 ) %>% # ONLY UP-EXPRESSED IN EXPERIMENTAL CNTRST (i.e. DOWN-EXP. IN Ref)
  group_by(Name) %>%
  summarise(across(sampleB, .fns = list), n = n()) 

QUERIES %>% dplyr::count(n)

QUERIES <- QUERIES %>% filter(n == 1) %>% unnest(sampleB)

QUERIES %>% dplyr::count(sampleB)

# ANNOT ====

orgdb <- "org.Hs.eg.db"

# semdata <- GOSemSim::godata(orgdb, ont="BP")
# write_rds(semdata, file = paste0(path, "hsGO_PB.rds"))

semdata <- read_rds(paste0(path, "hsGO_PB.rds"))

go_file <- paste0(path, 'Trinotate_GO_transcripts.txt')

MAP <- topGO::readMappings(go_file)

NAME2GO <- data.frame(Name = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERIES) %>%
  group_by(sampleB, Name) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())
  
gene2GO <- split(strsplit(NAME2GO$GO.ID, ";") , NAME2GO$Name)

gene2GO <- lapply(gene2GO, unlist)


print(which_sam <- QUERIES %>% distinct(sampleB) %>% pull())


allRes <- list()

for (i in which_sam) {
  
  cat("\nRunning GO enrichment for ", i)

  query.names <- NAME2GO %>% filter(sampleB %in% i) %>% distinct(Name) %>% pull()
  
  cat("\nUsing ", length(query.names), "transcrits")

  query.p <- RES.P %>% filter(Name %in% query.names) %>%
    group_by(Name) %>% sample_n(1) %>%
    pull(padj, name = Name)

  query.p <- query.p[match(query.names, names(query.p))]

  # identical(names(query.p), query.names)

  df <- GOenrichment(query.p, query.names, gene2GO, Nodes = 25, onto = "BP")

  allRes[[i]] <- data.frame(df, sampleB = i)
}

data <- do.call(rbind, allRes) %>% as_tibble()


