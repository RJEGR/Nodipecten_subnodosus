
library(tidyverse)

dir <- "~/Documents/GitHub/Nodipecten_subnodosus/data/"


f <- list.files(dir, pattern = "isoforms.counts.matrix", full.names = T) 

# cols <- read_tsv(f, skip = 1) %>% select(contains(".sorted.bam")) %>% names()

PCA <- function(f) {
  
  df <- read_tsv(f) 
  
  m <- df %>% select(contains(".isoforms.results")) %>% as("matrix")
  
  colNames <- colnames(m)
  
  colnames(m) <- gsub(".isoforms.results", "", colNames)
  
  data <- DESeq2::vst(round(m))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Method = basename(f))
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar),
    Method = basename(f))
  
  return(list(PCAdf, PCAvar))
}

PCAdf <- PCA(f)

percentVar <- PCAdf[[2]]$percentVar


PCAdf[[1]] %>%
  # mutate(Method = "Rnaspades") %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  mutate(g = substr(LIBRARY_ID, 1,2)) %>%
  ggplot(., aes(PC1, PC2, color = g)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  # scale_color_grey("") +
  facet_grid(~ Method) +
  # ylim(-100,100) + xlim(-100,100) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size = 7, alpha = 0.7) +
  geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 2.5, color = "black") +
  theme(legend.position = "top") -> p

p

# ggsave(p, filename = '03.Quantification.PCA.png', path = dir_out, width = 4, height = 3, device = png, dpi = 300)


PCAdf[[2]] %>%
  group_by(Method) %>%
  mutate(Method = "Rnaspades") %>%
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim), fill = Method, color = Method)) +
  geom_col(position = position_dodge2(), fill = "gray78", color = "gray78") +
  geom_line(aes(y = Varcum, group = Method)) +
  geom_point(aes(y = Varcum)) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "none", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 


# dendogram against first assembly

library(dendextend)

# Create two dendrograms
dend1 <- as.dendrogram (hc_regulatory_mir)
dend2 <- as.dendrogram (hc_pw_seq_align)

# Create a list to hold dendrograms
dend_list <- dendlist(dend1, dend2)


#  alignment quality vals
untangle(method = "ladderize") # 0.6520217
untangle(method = "step1side") #  0.3670336 <- most linear
untangle(method = "random", R = 10) # 0.4658022


# Align and plot two dendrograms side by side
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = FALSE,
    common_subtrees_color_lines = TRUE,
    common_subtrees_color_branches = TRUE) 


# Compute alignment quality. Lower value = good alignment quality
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>%
  entanglement() # Alignment quality

# Baker’s Gamma Index (see baker’s paper from 1974) is a measure of association (similarity) between two trees of Hierarchical clustering (dendrograms). It is defined as the rank correlation between the stages at which pairs of objects combine in each of the two tree