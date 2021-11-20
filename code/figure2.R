###
### FIGURE 2
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Scatterplot of spike only mean and variance
  ### x-axis is ln(mean)
  ### y-axis is ln(variance)
  ### 1 point for each spike (260 total)

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(GGally)
library(MASS)
library(reshape2)
library(viridis)
library(ggplot2)
library(zoo)
library(plyr)
library(dplyr)

#################
### ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### I/O
inputDir_v <- file.path("./data/spike_counts/25bp/ind/")
outDir_v <- file.path( "./output/fig2/")
meta_dt <- fread("./meta/meta.txt")

### Source functions
source(file.path("./code/functions.R"))

### Batch ### Samples S1 - S20 were used
batches_v <- "DNA160609LC" 

### Columns
cols_v <- c("V", "J", "spike.count")

#############
### INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Get input files
inputFiles_v <- list.files(inputDir_v)

### Name input files
inputNames_v <- paste0("S", gsub("^.*_S|.assemb.*", "", inputFiles_v))
names(inputFiles_v) <- inputNames_v

### Get desired samples from metadata
getSamples_v <- paste0("S", gsub("S", "", meta_dt[batch %in% batches_v, sample]))
inputFiles_v <- inputFiles_v[getSamples_v]

### Read in data
inputData_lsdt <- sapply(inputFiles_v, function(x) fread(file.path(inputDir_v, x), select = cols_v), simplify = F)

### Merge into one data.table
inputData_dt <- mergeDTs(inputData_lsdt, mergeCol_v = c("V", "J"))

####################
### CALCULATIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Get mean and variance
mean_v <- apply(inputData_dt[,mget(getSamples_v)], 1, mean)
var_v <- apply(inputData_dt[,mget(getSamples_v)], 1, var)

### Calculate dispersion
disp_lsdf <- apply(inputData_dt[,mget(getSamples_v)], 1, function(x) calcDisp(x, names = c("mu.nb", "theta")))
disp_df <- ldply(disp_lsdf, rbind)

### Combine
res_df <- data.frame(V = inputData_dt$V, J = inputData_dt$J, Mean = mean_v, Var = var_v)
res_df <- cbind(res_df, disp_df)

### Calculate logs
res_df$LogMean <- log(res_df$Mean)
res_df$LogVar <- log(res_df$Var)

### Calculate expected var
res_df$nb.Var.Estimate <- log(0.125) + 2*(res_df$mu.nb.Estimate)

############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

meanVar_gg <- ggplot(data = res_df, aes(x = LogMean, y = LogVar)) +
  geom_point() +
  geom_line(aes(x = mu.nb.Estimate, y = nb.Var.Estimate, col = "red")) +
  ggtitle("Mean vs. Variance") +
  labs(x = "ln Mean", y = "ln Variance") +
  guides(color = F) +
  big_label()

pdf(file = file.path(outDir_v, "fig2.pdf"))
print(meanVar_gg)
dev.off()