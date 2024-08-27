# red de datos de transcriptoma de coral

library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()


dir <- "~/Documents/GitHub/Nodipecten_subnodosus/Results/04.Quantification/"

count <- "CDHIT-95_good.Trinity.isoforms.counts.matrix$"

count_file <- list.files(path = dir, pattern = count, full.names = T)

count0 <- read.table(count_file, header=T, com='', row.names=1, check.names=F, sep='\t', stringsAsFactors = FALSE)

nr <- nrow(count0)

# Set out dir
out_dir <- paste0(dir, "wgcna_",Sys.Date())

system(paste0('mkdir -p ', out_dir))

setwd(out_dir)


# Filter data count ----
MIN_CPM=1;MIN_REPS=2

nrow(count0 <- count0[rowSums(edgeR::cpm(count0) > MIN_CPM) >= MIN_REPS, ])

nr

prevelancedf = apply(X = count0,
  MARGIN = 1,
  FUN = function(x){sum(x > 0)})


df = data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(count0)) %>% 
  as_tibble(rownames = 'id')

count <- round(count0)

LIBRARY_ID <-  gsub(".isoforms.results", "", basename(names(count))) 

LIBRARY_ID <- sapply(strsplit(LIBRARY_ID, "_CK"), `[`, 1)

names(count) <- LIBRARY_ID

conditions <- substr(LIBRARY_ID, 1,2)

conditions <- data.frame(conditions=factor(conditions), 
  LIBRARY_ID = factor(LIBRARY_ID))

rownames(conditions) <- colnames(count)

dds <- DESeq2::DESeqDataSetFromMatrix(count,
  conditions,
  design = ~LIBRARY_ID)

dds <- DESeq2::varianceStabilizingTransformation(dds)

datExpr <- SummarizedExperiment::assay(dds)

datExpr <- t(datExpr)

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

saveRDS(datExpr, file = paste0(dir, '/wgcna_datExpr.rds'))

# str(read_rds(paste0(dir, '/wgcna_datExpr.rds')))

# After detect the max power in the pickSoftThreshold.R ----

# The soft thresholding, is a value used to power the correlation of the genes to that threshold. The assumption on that by raising the correlation to a power will reduce the noise of the correlations in the adjacency matrix

# Blockwise construction ----

# Call the network topology analysis function

rds_f <- list.files(dir, pattern = "SoftThreshold", full.names = T)

sft <- readr::read_rds(rds_f)

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])
soft_values <- round(soft_values, digits = 2)
power_pct <- soft_values[which.max(soft_values)]
softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]
meanK <- sft$fitIndices[softPower,5]
softPower <- min(softPower)


allowWGCNAThreads()

# corOptionsList = list(use ='p', maxPOutliers= 0.05)

# if Test

str(datExpr <- read_rds(paste0(dir, '/wgcna_datExpr.rds')))

# dim(datExpr<- datExpr[1:36, sample(1:ncol(datExpr), size= 1000)])


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
