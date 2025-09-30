# Chord diagram, from Biological process to Sample groups


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "~/Documents/MANO_DELEON/DATOS_Paulina_dir/"

subdir <- ""

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

orgdb <- "org.Hs.eg.db"

semdata <- read_rds(paste0(dir, orgdb, ".PB.rds"))

go_file <- paste0(dir, 'Trinotate_transcripts.txt')

MAP <- topGO::readMappings(go_file)


RES.P <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds")) %>%
  filter(sampleA != "BLA_Reg" & !sampleB  %in% c("BLA_Reg", "LOL_Reg")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

paste_sam <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ':', collapse = ':')
}

QUERIES <- RES.P %>% 
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  distinct(ids, sampleX) %>%
  group_by(ids) %>%
  summarise(across(sampleX, .fns = paste_sam), n = n()) 

# Omit unigenes

QUERIES <- QUERIES %>% filter(n != 1) %>% unnest(sampleX)

QUERIES %>% dplyr::count(sampleX, sort = T) 


# First datavz

RES.P %>% 
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  distinct(ids, sampleX) %>%
  with(., table(ids, sampleX)) %>% data.frame() %>%
  ggplot(aes(y = ids, x = sampleX, fill = Freq)) +
  geom_tile()


my_cordplot <- function(DATA, grid.col, ...) {
  
  require(circlize)
  
  # dev.off()  
  
  circos.clear()
  
  DATA %>%
    chordDiagramFromDataFrame(
      directional = 1, 
      # col = colmat, 
      grid.col = grid.col,
      link.sort = T,
      # annotationTrack = "grid", 
      big.gap = 10, small.gap = 1,
      preAllocateTracks = list(track.height = 0.1),
      link.target.prop = FALSE)
}

dev.off()  

QUERIES %>%
  separate(sampleX, into = c("groupA", "groupB"), sep = ":") %>%
  # with(., table(groupA, groupB)) 
  sample_n(20) %>%
  my_cordplot()

dev.off()  

# second, using PB


NAME2GO <- data.frame(ids = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERIES) %>%
  group_by(sampleX, ids) %>%
  summarise(across(GO.ID, .fns = paste_go), n_gos = n()) %>%
  filter(n_gos > 1)

NAME2GO %>% dplyr::count(sampleX, sort = T) 
QUERIES %>% dplyr::count(sampleX, sort = T)

gene2GO <- split(strsplit(NAME2GO$GO.ID, ";") , NAME2GO$ids)

gene2GO <- lapply(gene2GO, unlist)


# Using Groups with GOs > 1,

which_sam <- NAME2GO %>% dplyr::count(sampleX, sort = T) %>%
 filter(n > 1) %>% distinct(sampleX) %>% pull()

str(which_sam)

# which_sam <- QUERIES %>% count(sampleX, sort = T) %>% filter(n > 1) %>% distinct(sampleX) %>% pull()

allRes <- list()

for (i in which_sam) {
  
  cat("\nRunning GO enrichment for ", i)
  
  query.names <- NAME2GO %>% filter(sampleX %in% i) %>% distinct(ids) %>% pull()
  
  cat("\nUsing ", length(query.names), "transcrits")
  
  query.p <- RES.P %>% filter(ids %in% query.names) %>%
    group_by(ids) %>% sample_n(1) %>%
    pull(padj, name = ids)
  
  query.p <- query.p[match(query.names, names(query.p))]
  
  # identical(names(query.p), query.names)
  
  df <- GOenrichment(query.p, query.names, gene2GO, Nodes = 25, onto = "BP")
  
  allRes[[i]] <- data.frame(df, sampleX = i)
}

data <- do.call(rbind, allRes) %>% as_tibble()

data %>% count(sampleX)


Chorddf <- data %>%
  distinct(Term, sampleX) %>%
  separate(sampleX, into = c("groupA", "groupB"), sep = ":", remove = F)

grid.nms <- c(unique(Chorddf$groupA),
  unique(Chorddf$groupB)
)

grid.col <- ggsci::pal_cosmic()(length(grid.nms))

grid.col <- structure(grid.col, names = sort(grid.nms))

grid.terms <- unique(Chorddf$Term)

grid.terms <- gsub(" ","\n",grid.terms)

grid.terms <- structure(rep("grey80", length(grid.terms)), names = grid.terms)

Chorddf %>% count(groupA, groupB)

dev.off() 
circos.clear()  

Chorddf %>%
  # with(., table(groupA, groupB))
  group_by(sampleX) %>% sample_n(2) %>%
  distinct(Term, groupA,groupB) %>%
  mutate(Term = gsub(" ","\n",Term)) %>%
  # my_cordplot(grid.col, )
  chordDiagramFromDataFrame(
    directional = -1, 
    grid.col = c(grid.col, grid.terms))

# Bind Unique GOs with goenrichment GOs


uniqueGenesSamples<- NAME2GO %>% dplyr::count(sampleX, sort = T) %>%
  filter(n == 1) %>% distinct(sampleX) %>% pull()

GOIDSdf <- NAME2GO %>%
  filter(sampleX %in% uniqueGenesSamples) %>%
  mutate(GO.ID = strsplit(GO.ID, ";")) %>%
  unnest(GO.ID) %>%
  distinct(GO.ID) %>%
  rbind(data %>% distinct(sampleX, GO.ID)) 

GOIDS <- GOIDSdf %>%
  # count(sampleX)
  pull(GO.ID, name = sampleX)

orgdb <- "org.Hs.eg.db"

semdata <- read_rds(paste0(dir, orgdb, ".PB.rds"))

REVIGO <- SEMANTIC_SEARCH(GOIDS, orgdb, semdata)

str(GOIDS)

dim(REVIGO)

dev.off() 
circos.clear()

REVIGO %>% 
  as_tibble(rownames = "GO.ID") %>%
  right_join(GOIDSdf) %>%
  distinct(parentTerm, sampleX) %>%
  separate(sampleX, into = c("groupA", "groupB"), sep = ":") %>%
  mutate(parentTerm = gsub(" ","\n",parentTerm)) %>%
  # with(., table(groupA, groupB))
  # group_by(sampleX) %>% sample_n(5) %>%
  my_cordplot()


# REVIGO %>% 
#   as_tibble(rownames = "GO.ID") %>% 
#   right_join(NAME2GO) %>% 
#   # distinct(parentTerm)
#   count(sampleX, parentTerm) %>%
#   drop_na() %>%
#   ggplot(aes(sampleX, parentTerm, fill = n)) +
#   geom_tile()



myf <- function(gene2GO) {
  
  
  # set <- names(gene2GO)
  
  golist <- gene2GO[[1]]
  
  cat("\n ", names(gene2GO), "\n ")
  
  cat("\n ", length(golist), "\n ")
  
  
  
  OUT <-  SEMANTIC_SEARCH(golist, orgdb, semdata)
  
  as_tibble(OUT) %>% mutate(Set = names(gene2GO))
}


myf(gene2GO[1])

OUT <- lapply(gene2GO, myf)


do.call(rbind, OUT) 

