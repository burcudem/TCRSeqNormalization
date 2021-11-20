###
### FIGURE 3
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
library(data.table)

#################
### ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### I/O

inputDir_v <- file.path("./data/spike_counts/25bp/ind/")
outDir_v <- file.path("./output/fig3/")
meta_dt <- fread(file.path("./meta/meta.txt"))

### Source functions
source(file.path("./code/functions.R"))

### Set divisions## samples S1 - S20 from DNA160609LC, samples 105-114 (new_dil) and 
#samples 125-134 (dil5) from were used from LIB170111LC ##

batches_v <- c("DNA160609LC", rep("LIB170111LC", 2))
dilutions_v <- c("original", "new_dil", "dil5")

### Assign reference and compare
refDil_v <- "original"
compDil_v <- c("new_dil", "dil5")

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

### List to hold data
data_lsdt <- list()

### Get each set of data
for (i in 1:length(batches_v)) {
  
  ## Get batch and dilution
  currBatch_v <- batches_v[i]
  currDil_v <- dilutions_v[i]
  name_v <- paste(currBatch_v, currDil_v, sep = "_")
  
  ## Get samples and files (Careful! This only works b/c there's no overlap in sample number between LIB170111LC and DNA160609LC...)
  currSamples_v <- paste0("S", meta_dt[batch == currBatch_v &
                                         dilution == currDil_v &
                                         is.na(tissue), sample])
  currFiles_v <- inputFiles_v[currSamples_v]
  
  ## Read in data and merge
  currData_lsdt <- sapply(currFiles_v, function(x) fread(file.path(inputDir_v, x), select = cols_v), simplify = F)
  currData_dt <- mergeDTs(currData_lsdt, mergeCol_v = c("V", "J"))
  
  ## Add to list
  data_lsdt[[name_v]] <- currData_dt
  
}

####################
### CALCULATIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### List to hold dispersion
disp_lsdf <- list()

### Calculate parameters for each set
for (i in 1:length(data_lsdt)) {
  
  ## Get data and info
  currName_v <- names(data_lsdt)[i]
  currBatch_v <- gsub("_.*$", "", currName_v)
  currDil_v <- gsub("^[A-Z0-9]+_", "", currName_v)
  currData_dt <- data_lsdt[[currName_v]]
  
  ## Calculate dispersion
  currSamples_v <- setdiff(colnames(currData_dt), c("V", "J"))
  currDispNames_v <- c(paste0("mu.nb_", currDil_v), paste0("theta_", currDil_v))
  currDisp_lsdf <- apply(currData_dt[,mget(currSamples_v)], 1, function(x) calcDisp(x, names = currDispNames_v))
  currDisp_df <- ldply(currDisp_lsdf, rbind)
  
  ## Add to list
  disp_lsdf[[currName_v]] <- currDisp_df
}

### Combine together and subset for desired columns
names(disp_lsdf) <- NULL
disp_df <- do.call(cbind, disp_lsdf)
getCols_v <- grep("Estimate|theta_", colnames(disp_df), value = T)
disp_df <- disp_df[,getCols_v]

### Add VJ column
vj_v <- paste(gsub("-", "", data_lsdt$DNA160609LC_original$V), data_lsdt$DNA160609LC_original$J, sep = "_")
disp_df$VJ <- vj_v


### This is not required because not doing the combined stuff, just using the 3 different sets
# ### Get difference between reference set (orig 20 spike only) and the compare sets
# refCol_v <- paste0('mu.nb_', refDil_v, ".Estimate")
# muDiff_lsv <- lapply(compDil_v, function(x) {
#   eCol_v <- paste0('mu.nb_', x, ".Estimate")
#   est <- ( sum(disp_df[[eCol_v]]) - sum(disp_df[[refCol_v]]) ) / 260
#   out <- disp_df[[eCol_v]] - est
#   return(out)
# })
# 
# ### Weight the differences and make new estimate
# disp_df$mu.nb_combined.Estimate <- (2*disp_df[[refCol_v]] + muDiff_lsv[[1]] + muDiff_lsv[[2]]) / 4

### Get ratios
disp_df$dil5.Ratio <- disp_df$mu.nb_dil5.Estimate - disp_df$mu.nb_original.Estimate
disp_df$new_dil.Ratio <- disp_df$mu.nb_new_dil.Estimate - disp_df$mu.nb_original.Estimate
old_mu_avg <- mean(disp_df$mu.nb_original.Estimate)
dil5_mu_avg <- mean(disp_df$mu.nb_dil5.Estimate)
new_dil_mu_avg <-mean(disp_df$mu.nb_new_dil.Estimate)
disp_df$old_avg_ratio <- log(disp_df$mu.nb_original.Estimate/old_mu_avg)
disp_df$dil5_avg_ratio <- log(disp_df$mu.nb_dil5.Estimate/dil5_mu_avg)
disp_df$new_dil_avg_ratio <- log(disp_df$mu.nb_new_dil.Estimate/new_dil_mu_avg)


### Subset for plot
#plot_df <- melt(disp_df[,c("VJ", "dil5.Ratio", "new_dil.Ratio")], id.vars = "VJ")
#MuDifferencesLabs <- c("1st set of 10 ST-only", "2nd set of 10 ST-only")
### Make boxplot for differences
#plot_gg <- ggplot(data = plot_df, aes(x = variable, y = value, fill = variable)) +
  #geom_boxplot() +
  #big_label() +
  #labs(y = expression(paste("Log Ratio")))+labs(x = NULL) +
  #guides(fill = F) + scale_x_discrete(labels=MuDifferencesLabs) + 
  #ggtitle(expression(paste("Log ratios of means: 10 ST relative to 20 ST")))
#pdf(file = file.path(outDir_v, "fig3.pdf"))
#print(plot_gg)
#dev.off()

##Correlation Plots##
cor1<-disp_df[,c("old_avg_ratio")]
cor2<-disp_df[,c("new_dil_avg_ratio")]
cor3<-disp_df[,c("dil5_avg_ratio")]
cormean<-disp_df[,c("old_avg_ratio", "new_dil_avg_ratio", "dil5_avg_ratio")]

# calculate correlations
meanCor1 <- cor.test(cor1,cor2)
meancorlm <- summary(lm(cor1 ~ cor2))

## Get results
meancor1_p <- meanCor1$p.value
meancor1_e <- meanCor1$estimate

## Get line paramters
meancorlmA_v <- meancorlm$coefficients[1,1]
meancorlmB_v <- meancorlm$coefficients[2,1]

## Make legend
meancor_Legend_v <- c(paste0("R = ", round(meancor1_e, digits = 2)))
#paste0("p = ", round(currP_v, digits = 1)))

plot(x = cormean$old_avg_ratio, xlab = "20 ST only",
     y = cormean$new_dil_avg_ratio, ylab = "2nd set of 10 ST only",
     main = "Scale Factors: 20 ST vs 2nd set of 10 ST", cex.lab=1, cex.main=1)
legend("bottomright", legend = meancor_Legend_v, bty = 'n', text.font=2, text.col = "Red" , c(as.expression(bquote(bold("Bold")))))
abline(a = meancorlmA_v, b = meancorlmB_v)

meanCor2 <- cor.test(cor1,cor3)
meancor2_e <- meanCor2$estimate
meancor_Legend_v2 <- c(paste0("R = ", round(meancor2_e, digits = 2)))
meancorlm2 <- summary(lm(cor1 ~ cor3))
meancorlmA_v2 <- meancorlm2$coefficients[1,1]
meancorlmB_v2 <- meancorlm2$coefficients[2,1]
 


plotcor1<-plot(x = cormean$old_avg_ratio, xlab = "20 ST only",
     y = cormean$dil5_avg_ratio, ylab = "1st set of 10 ST only",
     main = "Scale Factors: 20 ST vs 1st set of 10 ST", cex.lab=1, cex.main=1)
legend("bottomright", legend = meancor_Legend_v2, bty = 'n', text.font=2, text.col = "Red" , c(as.expression(bquote(bold("Bold")))))
abline(a = meancorlmA_v2, b = meancorlmB_v2)
print(plotcor1)

### Print
pdf(file = file.path(outDir_v, "fig3.pdf"), width = 9, height = 9)
print(plotcor1)
graphics.off()

cor_ratio1<-disp_df[,c("dil5.Ratio")]
cor_ratio2<-disp_df[,c("new_dil.Ratio")]
var(cor_ratio1)
var(cor_ratio2)
var.test(cor_ratio1, cor_ratio2)
var.test(cor1, cor2)
