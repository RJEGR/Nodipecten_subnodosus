# Run WGCNA by batch
# Ricardo Gomez-Reyes
# https://github.com/RJEGR/CORAL_PAPER

library(WGCNA)
library(DESeq2)
library(flashClust)
# library(DESeq2)

options(stringsAsFactors = FALSE)

MIN_REPS <- 2
MIN_CPM <- 1

dir <- "~/Documents/GitHub/Nodipecten_subnodosus/Results/04.Quantification/"

count_file <- list.files(path = dir, pattern = "CDHIT-95_good.Trinity.isoforms.counts.matrix", full.names = T)


nrow(count <- read.table(count_file, header=T, com='', 
  row.names=1, check.names=F, sep='\t', 
  stringsAsFactors = FALSE))

nrow(count <- count[rowSums(edgeR::cpm(count) > MIN_CPM) >= MIN_REPS, ])

LIBRARY_ID <-  gsub(".isoforms.results", "", basename(names(count))) 

LIBRARY_ID <- sapply(strsplit(LIBRARY_ID, "_CK"), `[`, 1)

names(count) <- LIBRARY_ID

conditions <- substr(LIBRARY_ID, 1,2)

conditions <- data.frame(conditions=factor(conditions), 
  LIBRARY_ID = factor(LIBRARY_ID))

rownames(conditions) <- colnames(count)

count <- round(count)

dds <- DESeqDataSetFromMatrix(count,
  conditions,
  design = ~LIBRARY_ID)

# dds <- DESeq(dds)

dds <- DESeq2::varianceStabilizingTransformation(dds)

datExpr <- assay(dds)

# barplot( colSums(datExpr))
# barplot( colSums(count))

plot(hclust(dist(t(datExpr))))
plot(hclust(dist(t(count))))

datExpr <- t(datExpr)

# datExpr <- t(head(datExpr, n = 100))

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK


if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}



max_power <- 30

powers = c(c(1:10), seq(from = 10, to = max_power, by=1)) 

allowWGCNAThreads()

cor_method =  "cor" # by default WGCNA::cor(method =  'pearson') is used, "bicor"

corOptionsList = list(use ='p') # maxPOutliers = 0.05, blocksize = 20000

sft <- pickSoftThreshold(datExpr, 
  powerVector = powers, 
  corFnc = cor_method,
  corOptions = corOptionsList,
  verbose = 5, 
  networkType = "unsigned")

rds_f <- paste0(dir, 'SoftThreshold_',cor_method, '.rds')

readr::write_rds(sft, file = rds_f)

sft <- readr::read_rds(rds_f)

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

power_pct <- soft_values[which.max(soft_values)]

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower,5]

hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

library(tidyverse)

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave

allowWGCNAThreads()

wd <- paste0(dir, ,"",Sys.Date())

system(paste0('mkdir ', wd))

setwd(wd)
# The variable datExpr now contains the expression data ready for network analysis.

setwd(dir)

bwnet <- blockwiseModules(datExpr, 
  maxBlockSize = 5000,
  power = softPower, 
  TOMType = "unsigned", 
  networkType = "unsigned",
  minModuleSize = 30,
  corType = "bicor",
  reassignThreshold = 0, 
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM-blockwise",
  verbose = 3)

saveRDS(bwnet, "bwnet.rds")

