# Ricardo Gomez-Reyes
# Visualize coverage of Nodipecten subnodosus

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)


dir <- "/Users/cigom/Documents/GitHub/Nodipecten_subnodosus/Results/Assembly-stats/Assembly-stats/cov_summaries"


f <- list.files(path = dir, pattern = "bowtie", full.names = T)

read_f <- function(f) { read_tsv(f) %>% mutate(Method = basename(f)) }

df1 <- lapply(f, read_f) # %>%  select(starts_with("PE")) %>% names()
df1 <- do.call(rbind, df1)

df1 <- df1 %>% mutate(Sample = sapply(strsplit(Sample, "_CK"), `[`, 1))

df1 %>% group_by(Method) %>% summarise(mean(overall_alignment_rate))


cols <- df1 %>%
  select(starts_with("paired_aligned")) %>% names()


df1 <- df1 %>% 
  pivot_longer(cols = all_of(cols), names_to = "Category", values_to =  "Reads") 

df1 <- df1 %>% mutate(Category = factor(Category, levels = cols))

df1 <- df1 %>%
  group_by(Sample, Method) %>%
  mutate(pct_align = Reads / sum(Reads))


category <- "paired_aligned_mate_none_halved"

p <- df1 %>%
  mutate(Method = stringr::str_to_title(Method)) %>%
  mutate(Method = gsub("_multiqc_bowtie2.Txt", "", Method)) %>%
  mutate(Align = ifelse(Category != category, "Transcriptome coverage", "Unalignment")) %>%
  # group_by(Sample, Method, Align) %>%
  # summarise(pct_align = sum(pct_align)) %>%
  filter(Align == "Transcriptome coverage") %>%
  ggplot(aes(y = Method, x = pct_align, fill = Category)) +
  geom_col() +
  facet_grid(  Sample ~., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) +
  scale_y_discrete(position = "right") +
  labs(y = "Assembly method", x = "% Alignment") +
  # coord_flip() +
  scale_fill_grey("") +
  # scale_fill_manual("Assembly method", values = c("black", "grey89")) +
  guides(fill=guide_legend(nrow = 10)) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  theme(legend.position = "right", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    axis.line.x = element_blank(),
    axis.line.y = element_blank())


ggsave(p, filename = 'alignment-methods.png', path = dir, width = 4.5, height = 8, device = png, dpi = 300)
