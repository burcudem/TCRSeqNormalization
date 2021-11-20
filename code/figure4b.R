###
### Figure 4B
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Boxplot/Stacked Scatter of monoclonal ratios before and after normalization
  ### Boxplot
    ### x-axis is norm/unnorm
    ### y-axis is monoclonal count
  ### Stacked scatter
    ### x-axis is sample ID
    ### y-axis is monoclonal count
    ### color is norm/unnorm

### Uses S73-S84 from LIB170111LC

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(data.table)
library(scales)
library(gridExtra)
library(grid)
library(reshape2)

#################
### ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### I/O

inputDir_v <- file.path("./data/normalized_clones/")
outDir_v <- file.path("./output/fig4/")
meta_dt <- fread(file.path("./meta/meta.txt"))

### Source functions
source(file.path("./code/functions.R"))

### Samples
batches_v <- "LIB170111LC"
samples_v <- paste0("S", 73:84)

### Columns
cols_v <- c("cloneId", "cloneCount", "cloneFraction", "nb.clone.count", "nb.clone.fraction", "V segments", "J segments", "aaSeqCDR3")

### Clones
ot1_lsv <- list("V" = "V121", "J" = "J2-7", "AA" = "CASSRANYEQYF")
p14_lsv <- list("V" = "V133", "J" = "J2-4", "AA" = "CASSDAGGRNTLYF")

#############
### INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Get input files and name
inputFiles_v <- list.files(inputDir_v)
names(inputFiles_v) <- gsub("^[A-Z0-9]*_|\\.clono.*", "", inputFiles_v)

### Subset
inputFiles_v <- inputFiles_v[names(inputFiles_v) %in% samples_v]

### Read in data
inputData_lsdt <- lapply(inputFiles_v, function(x){
  fread(file.path(inputDir_v, x), select = cols_v)
})

###########################
### CALCULATE DEVIATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################

### Empty variables
devData_dt <- NULL
dilutionData_dt <- data.table("sample" = paste0("S", 73:84), "dilution" = c(rep(1200, 3), rep(900, 3), rep(600, 3), rep(300, 3)))

### Calculate deviation for each sample
for (i in 1:length(inputData_lsdt)) {
  
  ##
  ## INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  ## Get name and data
  currSample_v <- names(inputData_lsdt)[i]
  currData_dt <- inputData_lsdt[[currSample_v]]
  
  ##
  ## GET CLONES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  ## Extract OT1
  currOT1_dt <- currData_dt[(aaSeqCDR3 == ot1_lsv$AA & 
                               `V segments` == ot1_lsv$V &
                               `J segments` == ot1_lsv$J),]
  currOT1_dt <- currOT1_dt[1,]
  
  ## Extract P14
  currP14_dt <- currData_dt[(aaSeqCDR3 == p14_lsv$AA &
                               `V segments` == p14_lsv$V &
                               `J segments` == p14_lsv$J),]
  currP14_dt <- currP14_dt[1,]
  
  ## Combine
  currMono_dt <- rbind(currOT1_dt, currP14_dt)
  
  # ## Add info
  # currMono_dt$Sample <- currSample_v
  # currMono_dt$Mono <- c("ot1", "p14")
  # currMono_dt$dilution <- dilutionData_dt[sample == currSample_v, dilution]
  
  ##
  ## CALCULATE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  ## Get deviation
  currDev_v <- apply(currMono_dt[,mget(c("cloneCount", "nb.clone.count"))], 2, function(x){
    abs((x[1] / (x[1] + x[2])) - 0.5)
  })
  
  ## Melt and add info
  currDev_dt <- convertDFT(melt(currDev_v), newName_v = "Count")
  currDev_dt$Sample <- currSample_v
  currDev_dt$dilution <- dilutionData_dt[sample == currSample_v, dilution]
  
  ## Add to overall table
  devData_dt <- rbind(devData_dt, currDev_dt)
  
} # for i

########################
### PREPARE FOR PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################

### Count versions
counts_v <- c("nb.clone.count", "cloneCount")

### Empty vars
dev_lsdt <- list()
sortDev_dt <- NULL

### Can't use normal melt() b/c need special ordering.
### Order based on nb.clone.count and then match others
for (i in 1:length(counts_v)) {
  
  ## Get name and data
  currName_v <- counts_v[i]
  currDev_dt <- devData_dt[Count == currName_v,]
  
  ## Sort
  if (currName_v == "nb.clone.count") {
    currDev_dt <- currDev_dt[order(value)]
  } else {
    currDev_dt <- currDev_dt[match(dev_lsdt$nb.clone.count$Sample, Sample),]
  } # fi
  
  ## Add index for ordering
  currDev_dt$index <- seq(i, (length(counts_v) * nrow(currDev_dt) - (length(counts_v)-i)), by = length(counts_v))
  
  ## Add to list
  dev_lsdt[[currName_v]] <- currDev_dt
  
  ## Add to sorted table
  sortDev_dt <- rbind(sortDev_dt, currDev_dt)
  
} # for i

############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Make stacked scatter
scatter_gg <- ggplot(data = sortDev_dt, aes(x = reorder(Sample, index), y = value, color = Count)) +
  geom_point() +
  ggtitle("Raw vs Normalized 50:50 Monoclonal \nP14 and OT-1 Samples") +
  labs(x = "Sample", y = "Deviation from Equal Ratio \nabs((A / A+B) - 0.5)") +
  big_label()

### Make box
box_gg <- ggplot(data = sortDev_dt, aes(x = Count, y = value, fill = Count)) +
  geom_boxplot() +
  ggtitle("Raw vs Normalized 50:50 Monoclonal \nP14 and OT1 Samples") + 
  guides(fill = F) +
  labs(x = NULL, y = "Deviation from Equal Ratio \nabs((A / A+B) - 0.5)") +
  #labels=c("Raw Clone Count", "Normalized Clone Count")+
  big_label()
box_gg

### Print
pdf(file = file.path(outDir_v, "fig4B_scatter.pdf"))
print(scatter_gg)
dev.off()

pdf(file = file.path(outDir_v, "fig4B_box.pdf"))
print(box_gg)
dev.off()