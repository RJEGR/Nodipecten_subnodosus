#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)

# if (length(argv) != 2) {
#   cat("Usage: rsem-plot-model sample_name output_plot_file\n")
#   q(status = 1)
# }

setwd("/Users/cigom/Documents/GitHub/Nodipecten_subnodosus/rsem_stats_dir")

argv <- list.files(path = getwd(), full.names = T, )

i <- 1

strvec <- strsplit(basename(argv[i]), split = "[.]")[[1]]

token <- strvec[1]

# stat.dir <- paste(argv[1], ".stat", sep = "")


stat.dir <- argv[i]

# if (!file.exists(stat.dir)) {
#   cat("Error: directory does not exist: ", stat.dir, "\n", sep = "")
#   q(status = 1)
# }

modelF <- paste(stat.dir, "/", token, ".model", sep = "")
cntF <- paste(stat.dir, "/", token, ".cnt", sep = "")

# pdf(argv[2])

con <- file(modelF, open = "r")	

# model type and forward probability
model_type <- as.numeric(readLines(con, n = 4)[1])  

# fragment length distribution
strvec <- readLines(con, n = 3)
vec <- as.numeric(strsplit(strvec[1], split = " ")[[1]])
maxL <- vec[2] # maxL used for Profile
x <- (vec[1] + 1) : vec[2]
y <- as.numeric(strsplit(strvec[2], split = " ")[[1]])
mode_len = which(y == max(y)) + vec[1]
mean <- weighted.mean(x, y)
std <- sqrt(weighted.mean((x - mean)^2, y))
plot(x, y, type = "h",
     main = "Fragment Length Distribution",
     sub = sprintf("Mode = %d, Mean = %.1f, and Std = %.1f", mode_len, mean, std),
     xlab = "Fragment Length",
     ylab = "Probability")
abline(v = mode_len, col = "red", lty = "dashed")

# mate length distribution
if (model_type == 0 || model_type == 1) bval <- as.numeric(readLines(con, n = 1)[1]) else bval <- 1

if (bval == 1) {
  list <- strsplit(readLines(con, n = 2), split = " ")
  vec <- as.numeric(list[[1]])
  maxL <- vec[2]
  x <- (vec[1] + 1) : vec[2]
  y <- as.numeric(list[[2]])
  mode_len = which(y == max(y)) + vec[1]
  mean <- weighted.mean(x, y)
  std <- sqrt(weighted.mean((x - mean)^2, y))
  plot(x, y, type = "h",
       main = "Read Length Distribution",
       sub = sprintf("Mode = %d, Mean = %.1f, and Std = %.1f", mode_len, mean, std),
       xlab = "Read Length",
       ylab = "Probability")
}
strvec <- readLines(con, n = 1)

# RSPD
bval <- as.numeric(readLines(con, n = 1)[1])
if (bval == 1) {
  bin_size <- as.numeric(readLines(con, n = 1)[1])
  y <- as.numeric(strsplit(readLines(con, n = 1), split = " ")[[1]])
  par(cex.axis = 0.7)
  barplot(y, space = 0, names.arg = 1:bin_size, main = "Read Start Position Distribution", xlab = "Bin #", ylab = "Probability")
  par(cex.axis = 1.0)
}
strvec <- readLines(con, n = 1)

# plot sequencing errors
if (model_type == 1 || model_type == 3) {
  # skip QD
  N <- as.numeric(readLines(con, n = 1)[1])
  readLines(con, n = N + 1)
  readLines(con, n = 1) # for the blank line
  
  # QProfile
  readLines(con, n = 1)

  x <- c()
  peA <- c() # probability of sequencing error given reference base is A
  peC <- c()
  peG <- c()
  peT <- c()
  
  for (i in 1 : N) {
    strvec <- readLines(con, n = 6)
    list <- strsplit(strvec[1:4], split = " ")

    vecA <- as.numeric(list[[1]])
    vecC <- as.numeric(list[[2]])
    vecG <- as.numeric(list[[3]])
    vecT <- as.numeric(list[[4]])

    if (sum(c(vecA, vecC, vecG, vecT)) < 1e-8) next
    x <- c(x, (i - 1))
    peA <- c(peA, ifelse(sum(vecA) < 1e-8, NA, -10 * log10(1.0 - vecA[1])))
    peC <- c(peC, ifelse(sum(vecC) < 1e-8, NA, -10 * log10(1.0 - vecC[2])))
    peG <- c(peG, ifelse(sum(vecG) < 1e-8, NA, -10 * log10(1.0 - vecG[3])))
    peT <- c(peT, ifelse(sum(vecT) < 1e-8, NA, -10 * log10(1.0 - vecT[4])))
  }

  matplot(x, cbind(peA, peC, peG, peT), type = "b", lty = 1:4, pch = 0:3, col = 1:4,
          main = "Observed Quality vs. Phred Quality Score",
          xlab = "Phred Quality Score",
          ylab = "Observed Quality")
  legend("topleft", c("A", "C", "G", "T"), lty = 1:4, pch = 0:3, col = 1:4)
} else {
  # Profile
  readLines(con, n = 1)

  x <- c()  
  peA <- c() # probability of sequencing error given reference base is A
  peC <- c()
  peG <- c()
  peT <- c()

  for (i in 1: maxL) {
    strvec <- readLines(con, n = 6)
    list <- strsplit(strvec[1:4], split = " ")

    vecA <- as.numeric(list[[1]])
    vecC <- as.numeric(list[[2]])
    vecG <- as.numeric(list[[3]])
    vecT <- as.numeric(list[[4]])

    if (sum(c(vecA, vecC, vecG, vecT)) < 1e-8) next
    x <- c(x, i)
    peA <- c(peA, ifelse(sum(vecA) < 1e-8, NA, (1.0 - vecA[1]) * 100))
    peC <- c(peC, ifelse(sum(vecC) < 1e-8, NA, (1.0 - vecC[2]) * 100))
    peG <- c(peG, ifelse(sum(vecG) < 1e-8, NA, (1.0 - vecG[3]) * 100))
    peT <- c(peT, ifelse(sum(vecT) < 1e-8, NA, (1.0 - vecT[4]) * 100))
  }

  matplot(x, cbind(peA, peC, peG, peT), type = "b", lty = 1:4, pch = 0:3, col = 1:4, main = "Position vs. Percentage Sequence Error", xlab = "Position", ylab = "Percentage of Sequencing Error")
  legend("topleft", c("A", "C", "G", "T"), lty = 1:4, pch = 0:3, col = 1:4)       
}

close(con)

# Alignment statistics
pair <- read.table(file = cntF, skip = 3, sep = "\t")

stat_len = dim(pair)[1]
upper_bound = pair[stat_len - 1, 1]
my_labels = append(0:upper_bound, pair[stat_len, 1])
my_heights = rep(0, upper_bound + 2)
dummy = sapply(1:(stat_len - 1), function(id) { my_heights[pair[id, 1] + 1] <<- pair[id, 2] })
my_heights[upper_bound + 2] = pair[stat_len, 2]
my_colors = c("green", "blue", rep("dimgrey", upper_bound - 1), "red")

barplot(my_heights, names.arg = my_labels,
        col = my_colors, border = NA,
        xlab = "Number of alignments per read",
        ylab = "Number of reads",
        main = "Alignment statistics")

pie_values = c(my_heights[1], my_heights[2], sum(my_heights[3:(upper_bound  + 1)]), my_heights[upper_bound + 2])
pie_names = c("Unalignable", "Unique", "Multi", "Filtered")
pie_labels = sprintf("%s %.0f%%", pie_names, pie_values * 100.0 / sum(pie_values))
par(fig = c(0.4, 1, 0.35, 0.95), new = T)
pie(pie_values, labels = pie_labels, col = c("green", "blue", "dimgrey", "red"), clockwise = T, init.angle = 270, cex = 0.8)
par(fig = c(0, 1, 0, 1))

dev.off.output <- dev.off()

# 

# Alignment statistics

read_cntF(argv[1])


read_cntF <- function(dir) {
  
  strvec <- strsplit(basename(dir), split = "[.]")[[1]]
  
  token <- strvec[1]

  stat.dir <- dir
  
  cntF <- paste(stat.dir, "/", token, ".cnt", sep = "")
  
  pair <- read.table(file = cntF, skip = 3, sep = "\t")
  
  stat_len = dim(pair)[1]
  
  upper_bound = pair[stat_len - 1, 1]
  my_labels = append(0:upper_bound, pair[stat_len, 1])
  my_heights = rep(0, upper_bound + 2)
  dummy = sapply(1:(stat_len - 1), function(id) { my_heights[pair[id, 1] + 1] <<- pair[id, 2] })
  my_heights[upper_bound + 2] = pair[stat_len, 2]
  
  pie_values = c(my_heights[1], my_heights[2], sum(my_heights[3:(upper_bound  + 1)]), my_heights[upper_bound + 2])
  pie_names = c("Unalignable", "Unique", "Multi", "Filtered")
  
  pie_frq <- pie_values * 100.0 / sum(pie_values)
  
  data.frame(pie_names, pie_values, pie_frq, "Sam" = token)
  
}

out <- lapply(argv, read_cntF)
out <- do.call( rbind, out)


out %>%
  ggplot2::ggplot(aes(x = Sam, y = pie_frq, fill = pie_names)) +
  ggplot2::geom_col() +
  theme_classic(base_size = 10, base_family = "GillSans") +
  theme(legend.position = "top")

out %>%
  ggplot2::ggplot(aes(x = Sam, y = pie_values, fill = pie_names)) +
  ggplot2::geom_col()


out %>% 
  select(-pie_frq ) %>%
  pivot_wider(names_from = pie_names, values_from = pie_values) %>% view()
