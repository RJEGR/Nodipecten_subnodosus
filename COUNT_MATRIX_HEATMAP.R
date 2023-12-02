#
# RUN SAMPLE-GROUP CLUSTERING 
# RICARDO GOMEZ-REYES
# 2023

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MANO_DELEON/"

# URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"
# 
# source(URL)

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

# NORMALIZE 

COUNTS <- round(COUNTS)

COUNTS <- DESeq2::vst(round(COUNTS))

# CLUSTERING SAMPLES
sample_cor = cor(COUNTS, method='pearson', use='pairwise.complete.obs')

# sample_dist = dist(sample_cor, method='euclidean')
# hc_samples = hclust(sample_dist, method='complete')
# hc_order <- hc_samples$labels[hc_samples$order]

hc <- heatmap(sample_cor, keep.dendro = T)

hc_samples <- as.hclust(hc$Rowv)

hc_order <- hc_samples$labels[hc$rowInd]

sample_cor %>% 
  as_tibble(rownames = 'Library_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(.colData) %>% 
  mutate(Library_ID = factor(Library_ID, levels = rev(hc_order))) -> sample_cor_long


# HEATMAP ====

library(ggh4x)

P <- sample_cor_long %>%
  ggplot(aes(x = Library_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggsci::scale_color_material(name = "", "blue-grey") +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1, size = 5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.05, "in"), label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(family = "GillSans", size = 7)))

P

# TOP PLOT ====

# recode_to <- c(`24` = "24 hpf", `110` = "110 hpf")

TOPDF <- sample_cor_long %>%
  distinct(Library_ID, Condition, Site) %>%
  # dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  # mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)

topplot <- TOPDF %>%
  ggplot(aes(y = y, x = Library_ID, color = Condition)) +
  geom_point(shape = 15, size = 2) +
  geom_text(aes(label = Site),  vjust = -1, hjust = 0.5, size = 1.5, family =  "GillSans", color = "black") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # see::scale_color_pizza(name = "", reverse = T) +
  # scale_color_manual("", values = c("#DADADA", "#D4DBC2")) + # "#4575b4", "#d73027"
  theme(legend.position = 'top',
    panel.border = element_blank(), 
    plot.background = element_rect(fill='transparent', color = 'transparent'),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank()) 
