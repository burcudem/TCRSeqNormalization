###
### FIGURE 4A
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Boxplots of spike distributions before and after normalization
  ### x-axis is sample ID
  ### y-axis is spike count (ln scale)

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
outDir_v <- file.path("./output/fig4/")
meta_dt <- fread(file.path("./meta/meta.txt"))
sf_dt <- fread(file.path("./meta/nb.scaling.factors.txt"))

### Source functions
source(file.path("./code/functions.R"))

### Samples
batches_v <- "DNA160609LC"

### Other
cols_v <- c("V", "J", "spike.count")
sfCols_v <- colnames(sf_dt)

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

### Remove '-' from V column
inputData_dt$V <- gsub("-", "", inputData_dt$V)

### Add duplicate rows for V121/V122
temp_dt <- inputData_dt[V == "V1212",]
temp_dt[, V := "V122"]
inputData_dt[V == "V1212", V := "V121"]
inputData_dt <- rbind(inputData_dt, temp_dt)

### Join V & J columns into one
inputData_dt$VJ <- paste(inputData_dt$V, inputData_dt$J, sep = "_")
inputData_dt <- inputData_dt[,mget(c("VJ", getSamples_v))]

sf_dt$VJ <- paste(sf_dt$V, sf_dt$J, sep = "_")
sf_dt <- sf_dt[,mget(c("VJ", setdiff(sfCols_v, c("V", "J"))))]

### Check
one_v <- setdiff(inputData_dt$VJ, sf_dt$VJ)
two_v <- setdiff(sf_dt$VJ, inputData_dt$VJ)
if (length(one_v) > 0 | length(two_v) > 0) stop("VJ combinations don't match between scaling factors and spike counts.")

### Convert integers to numeric
inputData_dt[, (getSamples_v) := lapply(.SD, as.numeric), .SDcols = getSamples_v]

#################
### NORMALIZE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### Create copy to normalize
normData_dt <- inputData_dt
normData_dt[1,1] <- normData_dt[1,1]

### Apply scaling factor to each (this is a pretty inefficient method...)
for (i in 1:nrow(sf_dt)) {
  
  ## Get data
  currVJ_v <- sf_dt[i,VJ]
  currSF_v <- sf_dt[i,scaling.factor]
  
  ## Apply for each sample
  for (j in 1:length(getSamples_v)) {
    
    ## Get sample
    currSample_v <- getSamples_v[j]
    
    ## Apply factor
    currNew_v <- normData_dt[VJ == currVJ_v, get(currSample_v)] * currSF_v
    currNew_v <- round(currNew_v, digits = 0)
    
    ## Add back
    normData_dt[VJ == currVJ_v, (currSample_v) := currNew_v]
    
  } # for j
  
} # for i

############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Melt
meltInput_dt <- melt(inputData_dt, id.vars = "VJ"); colnames(meltInput_dt) <- c("VJ", "Sample", "Raw")
meltNorm_dt <- melt(normData_dt, id.vars = "VJ"); colnames(meltNorm_dt) <- c("VJ", "Sample", "Norm")

### Merge
plot_dt <- merge(meltInput_dt, meltNorm_dt, by = c("VJ", "Sample"), sort = F)

### Melt
meltPlot_dt <- melt(plot_dt, id.vars = c("VJ", "Sample"))

### Apply ln
meltPlot_dt$ln <- log(meltPlot_dt$value)

### Make plot
plot_gg <- ggplot(data = meltPlot_dt, aes(x = Sample, y = ln, fill = variable)) +
  geom_boxplot() +
  labs(y = "ln(Count)", fill = "Count") +
  ggtitle("Raw and Normalized ST Count Distributions") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

### Print
pdf(file = file.path(outDir_v, "fig4A.pdf"), width = 9, height = 9)
print(plot_gg)
graphics.off()
