# Ricardo Gomez-Reyes
# Visualize Nx and BUSCO completness of Nodipecten subnodosus
# Calculate nx distribution


rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)


require(Biostrings)
require(dplyr)
require(ggplot2)

dir <- "/Users/cigom/Documents/GitHub/Nodipecten_subnodosus/Results/Assembly-stats/Assembly-stats/"


f <- list.files(dir, "fasta", full.names = T)

contig_Nx <- function(f) {
  
  DNA <- Biostrings::readDNAStringSet(f)
  
  contig_width <- sort(Biostrings::width(DNA), decreasing = T)
  
  # contig_df <- data.frame(width = contig_width, Assembly = basename(f))
  
  return(contig_width)
}


metrics_df <- function(f) {
  
  # width <- sort(seq(1,100), decreasing = TRUE)))
  
  width <- contig_Nx(f)
  
  assembly_seqs <- length(width)
  
  v <- seq(0.1,1, by = 0.1)
  
  Lx <- function(x) { sum(cumsum(width) < (sum(width) * x)) + 1 }
  
  l <- unlist( lapply(v, Lx))
  
  # add the number of sequences per Nx (field assembly_seqs)
  
  n_seqs <- sum(width > width[l[1]]) + 2
  
  # Cut the widths into chunks based on the Nx breakpoints 
  breakpts <- l
  
  # chunks <- cut(seq_along(width), breaks = l, labels = FALSE)
  
  chunks <- cut(rev(seq_along(width)), breaks = breakpts, labels = FALSE) 
  
  # Split the widths into chunks based on the cuts 
  chunks <- split(width, chunks) 
  
  n_seqs <- c(n_seqs,   as.vector(unlist(lapply(chunks, length))))
  
  n_frac <- n_seqs/assembly_seqs
  
  metrics_df <- data.frame(x = paste0("N", v*100), n = width[l], l = l, 
    n_seqs, n_frac, Assembly = basename(f))
  
  
  return(metrics_df)
  
}

# sum(metrics_df(f[1])$n_seqs)

df <- lapply(f, metrics_df)

df <- do.call(rbind, df)

df <- mutate(df, x = factor(x, levels = unique(df$x)))

# recode_to <- rev(unique(df$Assembly))

recode_to <- c("evigene_transcripts.fasta", "good.Trinity.fasta", "Trinity.fasta", "Trinity_ensamble_experimental.fasta", "Trinity_ensamble_bajas.fasta", 
  "bad.Trinity.fasta")


recode_to <- structure(c("Evigene", "good.Trinity","All-sam (Trinity)","Qphred_buenas (Pau)", "Qphred_bajas (Pau)", "bad.Trinity"), names = recode_to)

# recode_to <- structure(c("All-sam", "Qphred_buenas","Qphred_bajas","Evigene"), names = recode_to)

df <- mutate(df, Assembly = dplyr::recode_factor(Assembly, !!!recode_to))

# rnsps <- mean(contig_Nx(f[[1]]))
# trnt <- mean(contig_Nx(f[[2]]))

p <- ggplot(df, aes(x = x, y = n, group = Assembly, color = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  ggplot2::geom_path(linewidth = 1.5, lineend = "round") +
  geom_point(shape = 21, size = 4) +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  ggsci::scale_color_jco() +
  ggsci::scale_fill_jco() +
  # scale_color_grey("") +
  # scale_fill_grey("") +
  # scale_fill_manual("Assembly method", values = c("black", "grey89")) +
  guides(color=guide_legend(title = "", nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 
p
# annotate("text", y = rnsps, x = "N50", angle = 90, label = "label")

# p

ggsave(p, filename = 'Nx-methods.png', path = dir, width = 4, height = 3, device = png, dpi = 300)

ggplot(df, aes(y = x, x = n_frac, group = Assembly, fill = Assembly)) +
  # ggplot2::geom_col() +
  ggplot2::geom_col(position = position_dodge2(reverse = T)) +
  labs(x = "Frac. of Scaffolds", y = "Nx", fill = "Assembly method") +
  # scale_color_grey("") +
  # scale_fill_grey("") +
  ggsci::scale_color_jco() +
  ggsci::scale_fill_jco() +
  guides(fill=guide_legend(title = "", ncol = 1)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "left", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

df %>% group_by(Assembly) %>% summarise(n_seqs = sum(n_seqs), n_frac = sum(n_frac))
