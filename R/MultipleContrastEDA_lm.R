

# Evaluate what desing group is significant for the multiple contrast analysis
# Then re-doing deg analysis using 

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


DEGS_condition <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_condition.rds")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

DEGS <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds")) %>%
  filter(sampleA != "BLA_Reg" & !sampleB  %in% c("BLA_Reg", "LOL_Reg")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)


genes <- DEGS %>% distinct(ids) %>% pull(ids)

genes <- structure(genes, names = genes)

sum(keep <- rownames(datExpr) %in% sort(unique(names(genes))))

dim(datExpr <- datExpr[keep,])

genes <- genes[match(rownames(datExpr), names(genes))]

identical(rownames(datExpr), names(genes))

rownames(datExpr) <- genes


keep_cols <- colnames(datExpr) %in% levels(Manifest$Library_ID)

datExpr <- round(datExpr[,keep_cols])

DATA <- datExpr %>%  as_tibble(rownames = 'gene') %>% 
  pivot_longer(-gene, values_to = "expression", names_to = "Library_ID") %>%
  right_join(Manifest)
  # mutate(EDAD = as.integer(EDAD)) %>%
  # mutate_if(is.character, as.factor)



library(tidymodels)

# Wrap ------

fit_model <- function(data, formula) {
  
  
  require(tidymodels)
  
  # formula <- formula("expression ~ EDAD + gene")
  
  # formula <- c("expression ~ CONTRASTE_A + EDAD + gene")
  
  formula <- formula(formula)
  
  data_rec <- recipe(formula, data = data) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors())
  
  
  
  model_spec <- linear_reg() %>% set_engine("lm")
  
  data_workflow <- workflow() %>%
    add_recipe(data_rec) %>%
    add_model(model_spec)
  
  # Split training and test model set
  
  set.seed(123)
  
  data_split <- initial_split(data, prop = 0.75)
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  data_fit <- data_workflow %>% fit(data = data_train)
  
  # data_predictions <- data_fit %>%
  #   predict(new_data = data_test) %>%
  #   bind_cols(data_test)
  
  results <- data_fit %>%
    extract_fit_parsnip() %>%
    tidy()
  
  
  return(results)
}

# vars <- c("Condition","Site","design")

vars <- c("Condition","Site")

# fit_model(DATA, formula = "expression ~ Condition + Site")

fit_data <- list()

for(i in vars) {
  
  # form <- paste0("expression ~ ", i ," + gene")
  
  form <- paste0("expression ~ ", i)
  
  cat("\nCalculating ", form,"\n")
  
  fit_data[[i]] <- fit_model(DATA, form) %>% mutate(Intercept = i)
  
  # fit_data[[i]] <- OUT
}

fitdf <- do.call(rbind, fit_data)

fitdf <- fitdf %>%
  mutate(star = ifelse(p.value <.001, "***", 
    ifelse(p.value <.01, "**",
      ifelse(p.value <.05, "*", "ns"))))

fitdf

# recode_to <- c(  `Control` = "Control")

p <- fitdf %>%
  mutate(term = ifelse(Intercept == "Site" & grepl("Intercept", term), "Site_BLA", term)) %>%
  mutate(term = ifelse(Intercept == "Condition" & grepl("Intercept", term), "Condition_Bas", term)) %>%
  # dplyr::mutate(term = dplyr::recode_factor(term, !!!rev(recode_to))) %>%
  separate(term, into = c("facet", "term"), sep = "_") %>%
  mutate(x_star = estimate + (0.2+std.error)*sign(estimate)) %>%
  mutate(xmin = estimate-std.error, xmax = estimate+std.error) %>%
  ggplot(aes(y = term, x = estimate)) + # color = Intercept
  facet_grid(facet~ ., scales = "free", space = "free") +
  geom_point(size = 0.5,   color="gray40") +
  geom_text(aes(x = x_star, label = star),  
    vjust = 0.7, hjust = -0.3, size= 2, 
    # position=position_dodge(width = .5),
    family =  "GillSans") +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), 
    width = 0.1, alpha = 0.3, color = "gray40",
    # position=position_dodge(width = .5)
  ) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color = "black") +
  labs(
    y = "",
    x = "Effect Size") +
  scale_x_continuous(limits = c(-100,200)) +
  theme_bw(base_family = "GillSans", base_size = 12)  +
  theme(
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.text = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
    )

# p

ggsave(p, filename = 'MultipleContrast_effectSize.png', 
  path = dir, width = 2.5, height = 3, device = png, dpi = 700)



DEGS_condition

DataVizdf <- DEGS_condition %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2 ) %>%
  mutate(facet = ifelse( sign(log2FoldChange) == 1, "up in sampleA", "up in sampleB")) %>%
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  dplyr::count(sampleA, sampleB, facet, sampleX, sort = T) 

DataVizdf %>%
  # dplyr::mutate(facet = dplyr::recode_factor(sampleA, !!!recode_to)) %>%
  ggplot(aes(y = sampleA, x = sampleB, fill = n)) +
  facet_grid(~ facet, scales = "free") +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = n), size = 3, family = "GillSans", color = "white") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  labs(subtitle = "DEGS: FDR < 0.05 & abs(logFC) > 2") +
  theme(
    legend.position = "none",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 7),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 10)) -> P

P <- P +
  geom_segment(
    data = filter(DataVizdf, facet == "up in sampleA"), 
    x = 4, xend = 1, 
    y = 7, yend = 7, 
    colour = "gray7", 
    arrow = arrow(ends = "last", length = unit(0.15, "cm"))) +
  annotate("text", x = 4, y = 7, size = 3, label = "Contrast sence",  color = "gray7", family = "GillSans") +
  geom_segment(
    data = filter(DataVizdf, facet == "up in sampleB"), 
    x = 4, xend = 1, 
    y = 7, yend = 7, 
    colour = "gray7", 
    arrow = arrow(ends = "first", length = unit(0.15, "cm")))

P

ggsave(P, filename = 'MultipleContrast_heatmap.png', 
  path = dir, width = 7, height = 3, device = png, dpi = 700)


# Estimar La frecuencia de degs, en cada contraste, resumindo a los factores del disenio experimental
degs_df <- DataVizdf %>% group_by(sampleX) %>%  summarise(degs = sum(n), n_contrast = n()) %>% arrange(desc(degs))

DataVizdf %>%
  fit_model(formula = "n ~ sampleX") %>%
  filter(p.value < 0.05)

# Falta identificar si para cada grupo en sampleX, hay duplicados, 


UPSETDF <- DEGS_condition %>%
  # filter(padj < 0.05 & abs(log2FoldChange) > 2 ) %>%
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
  left_join(degs_df)

library(ggupset)

UPSETDF %>%
  # filter(n == 1) %>% ungroup() %>%
  ggplot(aes(x = sampleX, group = 1)) + # , fill = SIGN
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "degree") +
  theme_bw() +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans", base_size = 16) +
  # axis_combmatrix(levels = recode_to) +
  labs(x = '', y = 'Number of transcripts (up-expressed)') +
  # scale_color_manual("", values = col) +
  # scale_fill_manual("", values =  c("red", "black")) +
  guides(fill = guide_legend(title = "", nrow = 1)) -> p1

p1 <- p1 + theme(legend.position = "none",
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.background.y = element_blank())

p1


ggsave(p1, filename = 'UPSET.png', 
  path = dir, width = 7, height = 3, device = png, dpi = 700)


library(ggVennDiagram)

DF <- UPSETDF %>%  unnest(sampleX)

gene2ven <- split(DF$ids, DF$sampleX)

# gene2ven <- split(strsplit(UPSETDF$ids, "") , DF$Design)

gene2ven <- lapply(gene2ven, unlist)

str(gene2ven)

ggVennDiagram(gene2ven[1:7]) + scale_fill_gradient(low="grey90",high = "red")


