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

# TOPGO ====

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

write_rds(data, file = paste0(path, "DESEQ2TOPGO.rds"))

data <- read_rds(paste0(path, "DESEQ2TOPGO.rds"))

# REVIGO ====
# Additionally, run semantic simmilarity analysis in order to split parent terms themes

GO.IDS <- data %>% distinct(GO.ID) %>% pull() %>% sort()

SEMANTIC_SEARCH <- function(x, orgdb = "org.Ce.eg.db", semdata = semdata) {
  
  
  require(rrvgo)
  
  # semdata <- read_rds(paste0(wd, orgdb, ".rds"))
  
  x <- sort(x)
  
  SimMatrix <- calculateSimMatrix(x, 
    orgdb = orgdb,
    ont="BP", 
    semdata = semdata,
    method = 'Wang')
  
  data <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = orgdb) 
  
  y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
  
  data <- cbind(as.data.frame(y$points), data[match(rownames(y$points), data$go),])
  
  return(data)
}

# REVIGO <- SEMANTIC_SEARCH(GO.IDS, orgdb, semdata)

# data <- left_join(data, REVIGO, by = c("GO.ID" = "go"))

write_rds(data, file = paste0(path, "DESEQ2TOPGO.rds"))


# PLOT

data %>%
  group_by(sampleB) %>%
  mutate(size = size / max(size)) %>%
  filter(size > 0) %>%
  # arrange(desc(size), .by_group = T) %>%
  # mutate(Term = factor(Term)) %>%
  mutate(Term = fct_reorder(Term, size, .desc = F)) %>%
  # mutate(term = fct_reorder2(Term, sampleB, size, .desc = F)) %>%
  ggplot(aes(y = Term, x = size, color = -log10(p.adj.ks))) + # 
  facet_grid(sampleB ~ ., switch = "y", scales = "free_y") +
  geom_segment(aes(xend = 0, yend = Term), linewidth = 1) +
  labs(y = "Biological process (Up-expressed)", x = "Enrichment frac.") +
  # scale_color_viridis_c("-log10(padj)", option = "inferno") +
  theme(legend.position = "top",
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 90, size = 10)) -> p


ggsave(p, filename = 'GO_ENRICHMENT.png', path = path, width = 5, height = 10, device = png, dpi = 300)
