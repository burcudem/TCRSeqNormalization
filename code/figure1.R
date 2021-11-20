###
### FIGURE 1B
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Boxplots of spike distributions
  ### x-axis is spike identity
  ### y-axis is spike count/frequency
  ### 260 boxes that are made up of 10 or 20 points each from input samples

### Uses S1-S20 from DNA160609LC
  
####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(data.table)
library(ggplot2)

#################
### ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### I/O

inputDir_v <- file.path("./data/spike_counts/25bp/ind/")
outDir_v <- file.path(baseDir_v, "./output/fig1/")
meta_dt <- fread(file.path(baseDir_v, "./meta.txt"))

### Source functions
source(file.path("./code/functions.R"))

batches_v <- "DNA160609LC"

### Other
cols_v <- c("V", "J", "spike.count")

### Options - should y-axis be spike count or spike frequency (frequency means count / sum(count))
y_v <- "freq"
finalDir_v <- outDir_v

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

###############
### WRANGLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Remove '-' from V column and join V & J columns into one
inputData_dt$V <- gsub("-", "", inputData_dt$V)
inputData_dt$VJ <- paste(inputData_dt$V, inputData_dt$J, sep = "_")
inputData_dt <- inputData_dt[,mget(c("VJ", getSamples_v))]

### Convert to frequency, if needed
if (y_v == "freq") {
  for (col_v in getSamples_v) set(inputData_dt, j = col_v, value = inputData_dt[[col_v]] / sum(inputData_dt[[col_v]]))
  yAxis_v <- "ln(Spike Frequency)"
} else {
  yAxis_v <- "ln(Spike Count)"
}

### Find median of each spike
spikeMedian_v <- apply(inputData_dt[,mget(getSamples_v)], 1, median)
spikeMedian_dt <- data.table("VJ" = inputData_dt$VJ, "median" = spikeMedian_v)

### Sort by median
spikeMedian_dt <- spikeMedian_dt[order(median)]

### Melt spike count data
meltData_dt <- melt(inputData_dt, id.vars = "VJ")

### Turn VJ into a factor with levels taken from median table
meltData_dt$VJ <- factor(meltData_dt$VJ, levels = spikeMedian_dt$VJ)

### Take log of spike count
meltData_dt$ln <- log(meltData_dt$value)

#############
### PLOTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Raw
pdf(file = file.path(finalDir_v, "fig1_raw.pdf"), width = 9, height = 9)
ggplot(data = meltData_dt, aes(x = VJ, y = value)) +
  geom_boxplot() +
  my_theme() + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank()) +
  labs(x = "Spike", y = "Spike Frequency") +
  ggtitle("Distribution and Variation of Spike Count Frequency")
graphics.off()

### Nat. Log
pdf(file = file.path(finalDir_v, "fig1_ln.pdf"), width = 9, height = 9)
#lty<-c("y=ln(1/260)")
p<-ggplot(data = meltData_dt, aes(x = VJ, y = ln)) +
  geom_boxplot() +
  my_theme() + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank()) + 
  labs(x = "ST", y = "ln(ST Frequency)") +
  ggtitle("Distribution and Variation of ST Count Frequency") 
p + geom_hline(aes(yintercept=(log(1/260)), col="red"))+geom_text(aes(0,(log(1/260)),label = "                       ln(1/260)", vjust=-1)) + theme(legend.position = "none")
graphics.off() 


####################
### DATA SUMMARY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

###############
### PREPARE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Get vectors to iterate over
allST_v <- unique(meltData_dt$VJ)
allSamples_v <- unique(meltData_dt$variable)

### Make empty vectors to fill
rawSampleVar_v <- lnSampleVar_v <- numeric(length(allST_v))
rawStVar_v <- lnStVar_v <- numeric(length(allSamples_v))

rawSampleIQR_v <- lnSampleIQR_v <- numeric(length(allST_v))
rawStIQR_v <- lnStIQR_v <- numeric(length(allSamples_v))

##################################
### SAMPLE TO SAMPLE VARIATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################

### Calculate variation between 20 samples for 260 ST
for (i in 1:length(allST_v)) {
  
  ## Get ST and data
  currST_v <- allST_v[i]
  currData_dt <- meltData_dt[VJ == currST_v,]
  
  ## Check for infinite
  isInfinite_v <- which(is.infinite(currData_dt$ln))
  if (length(isInfinite_v) > 0) {
    infSamples_v <- as.character(currData_dt[isInfinite_v,variable])
    cat(sprintf("ST %s has Inf value for sample(s): %s", currST_v, paste(infSamples_v, collapse = "; ")))
    for (row_v in isInfinite_v) currData_dt[row_v, ln := NA]
  }
  
  ## Calculate SD
  currRawSD_v <- sd(currData_dt$value, na.rm = T)
  currLnSD_v <- sd(currData_dt$ln, na.rm = T)
  ## Calculate IQR
  currRawIQR_v <- IQR(currData_dt$value, na.rm = T)
  currLnIQR_v <- IQR(currData_dt$ln, na.rm = T)
  
  ## Add to vector
  rawSampleVar_v[i] <- currRawSD_v
  lnSampleVar_v[i] <- currLnSD_v
  
  rawSampleIQR_v[i]=currRawIQR_v
  lnSampleIQR_v[i]=currLnIQR_v
  
} # for i

### Get mean SD and sd of SD
rawMeanSampleVar_v <- mean(rawSampleVar_v)
rawSDSampleVar_v <- sd(rawSampleVar_v)

lnMedianIQR_v<-median(lnSampleIQR_v)
rawMedianIQR_v<-median(rawSampleIQR_v)

print(lnMedianIQR_v)
print(rawMedianIQR_v)

lnMeanSampleVar_v <- mean(lnSampleVar_v)
lnSDSampleVar_v <- sd(lnSampleVar_v)

##########################
### ST TO ST VARIATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################

### Calculate variation between 260 spikes for 20 samples
for (i in 1:length(allSamples_v)) {
  
  ## Get sample and data
  currSample_v <- allSamples_v[i]
  currData_dt <- meltData_dt[variable == currSample_v, ]
  
  ## Check for infinite
  isInfinite_v <- which(is.infinite(currData_dt$ln))
  if (length(isInfinite_v) > 0) {
    infVJ_v <- as.character(currData_dt[isInfinite_v, VJ])
    cat(sprintf("Sample %s has Inf value for ST(s): %s", currSample_v, paste(infVJ_v, collapse = "; ")))
    for (row_v in isInfinite_v) currData_dt[row_v, ln := NA]
  } # fi
  
  ## Calculate SD
  currRawSD_v <- sd(currData_dt$value, na.rm = T)
  currLnSD_v <- sd(currData_dt$ln, na.rm = T)
  currLnIQR_v <- IQR(currData_dt$ln, na.rm = T)
  currRawIQR_v <- IQR(currData_dt$value, na.rm = T)
  
  
  ## Add to vector
  rawStVar_v[i] <- currRawSD_v
  lnStVar_v[i] <- currLnSD_v
  lnStIQR_v[i]=currLnIQR_v
  rawStIQR_v[i]=currRawIQR_v
  
} # for i

### Get mean SD and sd of SD
rawMeanSTVar_v <- mean(rawStVar_v)
rawSDSTVar_v <- sd(rawStVar_v)

lnMedianIQR_v<-median(lnSampleIQR_v)

lnMeanSTVar_v <- mean(lnStVar_v)
lnSDSTVar_v <- sd(lnStVar_v)

lnMedianSTIQR_v<-median(lnStIQR_v)
rawMedianSTIQR_v<-median(rawStIQR_v)

print(lnMedianSTIQR_v)
print(rawMedianSTIQR_v)
###############
### COMBINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Make table
out_mat <- matrix(nrow = 2, ncol = 6)

### Add sample to sample
out_mat[1,] <- c(rawMeanSampleVar_v, rawSDSampleVar_v,
                 lnMeanSampleVar_v, lnSDSampleVar_v, 
                 rawMedianIQR_v, lnMedianIQR_v)

### Add ST to ST
out_mat[2,] <- c(rawMeanSTVar_v, rawSDSTVar_v,
                 lnMeanSTVar_v, lnSDSTVar_v,
                 rawMedianSTIQR_v,lnMedianSTIQR_v)

### Convert to data.table
rownames(out_mat) <- c("Sample_to_Sample", "ST_to_ST")
colnames(out_mat) <- c("Mean Raw SD", "SD of Raw SD", "Mean LN SD", "SD of LN SD", "Median Raw IQR", "Median LN IQR")
out_dt <- convertDFT(out_mat, newName_v = "Type")

### Write output
write.table(out_dt,
            file = file.path(outDir_v, "fig1_data.txt"),
            sep = '\t', quote = F, row.names = F)
