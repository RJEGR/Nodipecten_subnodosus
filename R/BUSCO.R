
rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)


require(tidyverse)


dir <- "/Users/cigom/Documents/GitHub/Nodipecten_subnodosus/Results/Assembly-stats/Assembly-stats/"

f2 <- list.files(dir, pattern = "Assembly-BUSCO.tsv", full.names = T)


df2 <- read_tsv(f2) %>%
  mutate(my_species = gsub("_odb10", "", my_species)) %>%
  mutate(my_species = stringr::str_to_title(my_species)) %>%
  mutate(Method = stringr::str_to_title(Method))

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

#names <- c("Complete", "Duplicated", "Fragmented", "Missing")

names <- c("S", "D", "F", "M")

labels <- c("Complete (C) and single-copy (S)",
  "Complete (C) and duplicated (D)",
  "Fragmented (F)  ",
  "Missing (M)")

col <- structure(col, names = rev(names))

my_sp_lev <- c("Mollusca","Metazoa","Eukaryota","Mammalia","Bacteria")

figure <- df2 %>%
  mutate(facet = "Completeness") %>%
  filter(Method != "Trinity-Longest") %>%
  # filter(my_species != "Bacteria") %>%
  mutate(category = factor(category, levels = rev(names))) %>%
  mutate(my_species = factor(my_species, levels = my_sp_lev)) %>%
  ggplot(aes(x = Method, y = my_percentage, fill = category)) +
  facet_grid(my_species~ facet) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "Assembly method") +
  coord_flip() +
  scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

ggsave(figure, filename = 'BUSCO-methods.png', path = dir, width = 5, height = 6, device = png, dpi = 300)
