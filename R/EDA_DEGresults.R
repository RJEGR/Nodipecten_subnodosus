
# Ricardo Gomez-Reyes
# Process DE results and identify intersection of experimental groups using venn diagrams
# According to degs filter count matrix and calculate zscore to previs patterns expressiosn per experimental groups)
# PLOT PCA
# PLOT Heatmap 

# install.packages("ggVennDiagram")


rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)
library(ggVennDiagram)

path <- "~/Documents/MANO_DELEON/DATOS_Paulina_dir/"

Manifest <- "metadata_experimental.tsv"

Manifest <- read_tsv(list.files(path = path, pattern = Manifest, full.names = T)) %>%
  dplyr::rename("LIBRARY_ID" = "Library_ID") %>%
  mutate(groups = paste(Site, Condition, sep = "_"))


f <- list.files(path = path, pattern = "DEG_", full.names = T)

PCA <- function(datExpr) {
  
  # require(DESeq2)
  
  data <- DESeq2::varianceStabilizingTransformation(round(datExpr))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar)
    # Method = basename(f)
  )
  
  return(list(PCAdf, PCAvar))
}

read_csv <- function(f) {
  
  out <- data.frame(fileName = basename(f), read.csv(f)) %>% as_tibble()
  return(out)
}

DF <- lapply(f, read_csv)

DF <- do.call(rbind, DF)

unique(DF$fileName)

# Aclarar grupos con paulina



recode_ref <- c("DEG_Bas_LOL_Ref.csv" = "LOL",
  "DEG_Cao_Bas_BLA.csv" = "Basal",
  "DEG_Cao_Bas_LOL.csv" = "Basal",
  "DEG_Cao_LOL_Ref.csv" = "LOL",
  "DEG_Gral_Cao_Bas_Ref.csv" = "Basal",
  "DEG_Cte_Bas_LOL.csv" = "Basal",
  "DEG_Cte_Bas_BLA.csv" = "Basal",
  "DEG_Cte_LOL_Ref" = "LOL",
  "DEG_Gral_Cte_Bas_Ref.csv" = "Basal")

recode_contrast <- c("DEG_Bas_LOL_Ref.csv" = "BLA",
  "DEG_Cao_Bas_BLA.csv" = "Chaotic",
  "DEG_Cao_Bas_LOL.csv" = "Chaotic",
  "DEG_Cao_LOL_Ref.csv" = "BLA",
  "DEG_Gral_Cao_Bas_Ref.csv" = "Chaotic",
  "DEG_Cte_Bas_LOL.csv" = "Cte",
  "DEG_Cte_Bas_BLA.csv" = "Cte",
  "DEG_Cte_LOL_Ref" = "BLA",
  "DEG_Gral_Cte_Bas_Ref.csv" = "Cte")



recode_design <- c("DEG_Bas_LOL_Ref.csv" = "Differences in populations (Basal site)",
  "DEG_Cao_Bas_BLA.csv" = "Specific acclim. to chaotic challenge (BLA)",
  "DEG_Cao_Bas_LOL.csv" = "Specific acclim. to chaotic challenge (LOL)",
  "DEG_Cao_LOL_Ref.csv" = "Differences in populations (Chaotic treatment)",
  "DEG_Gral_Cao_Bas_Ref.csv" = "Global acclimatation response to challenge (Chaotic)",
  "DEG_Cte_Bas_LOL.csv" = "Specific acclim. to cte challenge (LOL)",
  "DEG_Cte_Bas_BLA.csv" = "Specific acclim. to cte challenge (BLA)",
  "DEG_Cte_LOL_Ref" = "Differences in populations (Cte treatment)",
  "DEG_Gral_Cte_Bas_Ref.csv" = "Global acclimatation response to challenge (Cte)")

DF <- DF %>% dplyr::mutate(sampleA = dplyr::recode(fileName, !!!recode_ref)) 
DF <- DF %>% dplyr::mutate(sampleB = dplyr::recode(fileName, !!!recode_contrast)) 
DF <- DF %>% dplyr::mutate(Design = dplyr::recode(fileName, !!!recode_design)) 


DF <- DF %>% 
  mutate(sampleX = ifelse(sign(log2FoldChange) == 1, sampleA, sampleB)) 

DF %>% 
  count(fileName, Design, sampleA, sampleB, sampleX) %>% view()

DF <- DF %>%
  filter(fileName %in% c("DEG_Cao_Bas_BLA.csv", "DEG_Cao_Bas_LOL.csv", "DEG_Cte_Bas_LOL.csv", "DEG_Cte_Bas_BLA.csv")) %>%
  filter(sampleX != "Basal")

DF %>% 
  count(sampleB) 
# 
# 

gene2ven <- split(strsplit(DF$transcript, ";") , DF$Design)

gene2ven <- lapply(gene2ven, unlist)

ggVennDiagram(gene2ven) + scale_fill_gradient(low="grey90",high = "red")

ggVennDiagram(gene2ven[1:2])
ggVennDiagram(gene2ven[3:4])

# query_groups <- c("BLA", "LOL")
  
# q <- which(names(gene2ven) %in% query_groups)

# ggVennDiagram(gene2ven[q])


# step 2 -----
# PCA -----
str(query <- unique(DF$transcript))

f <- list.files(path = path, pattern = "isoforms.counts_experimental.matrix", full.names = T)

dim(M <- read.delim(f, row.names = 1))
keep <-  rownames(M) %in% query
dim(M <- M[keep,])

keep <-  !grepl("Reg", colnames(M))
dim( M <- M[,keep])

M <- as(M, "matrix")

PCAdf <- PCA(M)

percentVar1 <- paste0(PCAdf[[2]]$percentVar[1], "% ,",PCAdf[[1]]$percentVar[1])
percentVar1 <- paste0("PC1, VarExp: ", percentVar1, "%")

percentVar2 <- paste0(PCAdf[[2]]$percentVar[2], "% ,",PCAdf[[1]]$percentVar[2])
percentVar2 <- paste0("PC2, VarExp: ", percentVar2, "%")


PCAdf_ <- PCAdf[[1]]


PCAdf_ %>%
  mutate(LIBRARY_ID = rownames(.)) %>% 
  left_join(Manifest, by = "LIBRARY_ID") %>%
  ggplot(., aes(PC1, PC2)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  xlab(percentVar1) + ylab(percentVar2) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = HPF), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme(
    panel.grid = element_blank(),
    legend.position = 'top',
    legend.key.width = unit(0.2, "mm"),
    legend.key.height = unit(0.2, "mm")) -> p

p + 
  geom_point(aes(color = Condition), shape = 1, size = 7) + 
  ggrepel::geom_text_repel(aes(label = Site), family = "GillSans")



# Heatmap w/ zscores -----

z_scores <- function(x) {(x-mean(x))/sd(x)}

# barplot(colSums(M))

M <- DESeq2::varianceStabilizingTransformation(round(M))


# barplot(colSums(M))

# If summarise

groups <- paste(Manifest$Site, Manifest$Condition, sep = "_")

groups <- groups[match(Manifest$LIBRARY_ID, colnames(M))]

identical(Manifest$LIBRARY_ID, colnames(M))

row_group_means <- function(M, groups) {
  sapply(unique(groups), function(g) {
    rowMeans(M[, groups == g, drop=FALSE])
  })
}


# M <- row_group_means(M, groups)

# Cont. .
M <- apply(M, 2, z_scores)

heatmap(M)

h <- heatmap(M, keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
# plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]


# DATA <- data.frame(M) %>% 
#   as_tibble(rownames = 'LIBRARY_ID') %>%
#   pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "gene_id") %>%
#   left_join(Manifest, by = "LIBRARY_ID") %>%
#   mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))

names_to <- "LIBRARY_ID" # "LIBRARY_ID"

DATA <- data.frame(M) %>% 
  as_tibble(rownames = 'gene_id') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = names_to) %>%
  left_join(Manifest, by = names_to) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))
  # mutate(groups = factor(groups, levels = rev(hc_order))) %>%
  # left_join(distinct(Manifest, Condition, Site, groups), by = names_to) 

# If additional labels are provided, such as uniprot, etc.

# labels <- .TARGETDB %>% select(gene_id, description)
# arrange(desc(Name)) %>%
# group_by(Name) %>%
# mutate(Label = Name, row_number = row_number(Label)) %>%
# mutate(Label = factor(paste(Label, row_number, sep = "__"), 
#   levels = rev(paste(Label, row_number, sep = "__"))))


# labels <- structure(labels$description, names = labels$gene_id)

labels <- structure(rownames(M), names = rownames(M)) 

labels <- labels[match(order_genes, names(labels))]

identical(order_genes, names(labels))

order_genes <- labels

order_genes <- gsub("TRINITY_*", "", order_genes)

lo = floor(min(DATA$fill))
up = ceiling(max(DATA$fill))
mid = (lo + up)/2


library(ggh4x)

DATA %>%
  ggplot(aes(x = groups, y = gene_id, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.025, width = 1) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  # guides(y.sec = guide_axis_manual(labels = order_genes, label_size = 5, label_family = "GillSans")) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> P


# P

P <- P + 
  ggh4x::facet_nested(~ Condition+Site, nest_line = T, scales = "free_x", space = "free_x") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = 'white', color = 'white')) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(2.5, "in"),
    barheight = unit(0.07, "in"), 
    label.position = "top",
    title = "Row z-score",
    title.position  = "top",
    title.theme = element_text(size = 7, family = "GillSans", hjust = 1),
    alignd = 0.8,
    ticks.colour = "black", ticks.linewidth = 0.25,
    frame.colour = "black", frame.linewidth = 0.25,
    label.theme = element_text(size = 7, family = "GillSans")))

P

dir_out <-"~/Documents/MANO_DELEON/DATOS_Paulina_dir/PUBLICATION_DIR"

ggsave(P, filename = 'Zscore_hetmap.png', path = dir_out, 
  width = 8, height = 8, device = png, dpi = 300)
