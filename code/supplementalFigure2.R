library(dplyr)
library(data.table)
library(magrittr)
library(devtools)
#library(plot3D)
library(latticeExtra)

### Make a 3D Barplot of normalized clone counts from mixcr output files.

# This function takes 4 parameters:
# clone.dir is the path to a directory containing exported clones from MiXCR, with normalized counts
# out.dir is path of desired output directory for formatted tables and vj removal comparison
# spike.file is path to text_barcodesvj.txt
# metadata.file is a file containing desired clone files in 1st column, treatment in second column
# batch is name of batch to subset from metadata. Requires that metadata.file has column called "batch"

### Arguments
clone.dir <- file.path("./data/normalized_clones/")
out.dir <- file.path("./output/suppFig2/")
spike.file <- file.path("./meta/text_barcodesvj.txt")
metadata.file <- file.path("./meta/meta.txt")

### Source functions
source(file.path("code/functions.R"))

###
### Spike Count Table ~~~~~~~~~~~~~~~~~~~~~~~~~
###

# Read in spike count table
spike.table <- fread(spike.file)

# Remove trailing dashes for hygiene's sake. Also change V1212 (for consistency)
spike.table$V <- gsub("-","", spike.table$V)

# Extract the V and J sequences from spike table
vsegs <- unique(spike.table$V)
jsegs <- unique(spike.table$J)

# Read in metadata
meta_dt <- fread(metadata.file)

###
### Clone Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

# Read in normalized clone files and order them
all.clones <- list.files(clone.dir)
all.clones <- all.clones[order(as.numeric(gsub(".*_S|_alignment_.*|\\.clonotypes.*", '', all.clones)))]

### Columns to subset for
segments_v <- c("V segments", "J segments")
freqCol_v <- c("cloneFraction", "nb.clone.fraction")

###
### Create matrix for 3D plotting ~~~~~~~~~~~~~~~
###

for (j in 1:length(freqCol_v)) {
  
  ## Get columns
  currFreqCol_v <- freqCol_v[j]
  currCols_v <- c(segments_v, currFreqCol_v)
  cat(sprintf("Working on %s\n", currFreqCol_v))
  
  ## Make directory
  currOut_v <- mkdir(out.dir, currFreqCol_v)
  
  for (i in 1:length(all.clones)) {
    
    ## Get file and data
    currFile_v <- all.clones[i]
    currData_dt <- fread(file.path(clone.dir, currFile_v))
    
    ## Get names
    currBatch_v <- strsplit(currFile_v, split = "_|\\.")[[1]][1]
    currSample_v <- strsplit(currFile_v, split = "_|\\.")[[1]][2]
    
    ## Get monoclonal
    currMono_v <- meta_dt[batch == currBatch_v & sample == gsub("^S", "", currSample_v), monoclonal]
    
    ## Make output name
    if (is.na(currMono_v)) {
      currName_v <- paste(currBatch_v, currSample_v, "3d.pdf", sep = "_")
    } else {
      currName_v <- paste(currBatch_v, currSample_v, currMono_v, "3d.pdf", sep = "_")
    } # fi
    
    ## Remove V and J that we don't have
    currSub_dt <- currData_dt[`V segments` %in% vsegs,]
    currSub_dt <- currSub_dt[`J segments` %in% jsegs,]
    
    ## Subset columns
    currSub_dt <- currSub_dt[,mget(currCols_v)]
    
    ## Make count table
    currVJCount_dt <- currSub_dt[, .("count" = sum(get(currFreqCol_v))), by = segments_v]
    colnames(currVJCount_dt) <- c("V", "J", "count")
    
    ## Full version to hold 0's
    currFullVJCount_dt <- spike.table[,mget(c("V", "J"))]
    currFullVJCount_dt$count <- 0
    
    ## Combine, remove extra column, convert NA to 0
    currPlot_dt <- merge(currFullVJCount_dt, currVJCount_dt, by = c("V", "J"), all.x = T, sort = F)
    currPlot_dt <- currPlot_dt[, count.x := NULL]
    colnames(currPlot_dt) <- c("V", "J", "count")
    currPlot_dt[is.na(count), count := 0]
    
    ## determine if should take the log
    if (length(grep("[Cc]ount", currFreqCol_v)) == 0) {
      y <- currPlot_dt$count
      currLab_v <- currFreqCol_v
    } else {
      y <- log2(currPlot_dt$count)
      currLab_v <- paste0("log2(", currFreqCol_v, ")")
    }
    
    ## Assign for easier use
    #y <- currPlot_dt$count
    x <- factor(currPlot_dt$J, levels = unique(currPlot_dt$J))
    z <- factor(currPlot_dt$V, levels = unique(currPlot_dt$V))
    d <- data.frame(x, y, z)
    
    ## Set colors
    my.colors <- rep(c("skyblue", "blue4", "darkolivegreen1", "darkolivegreen4","darksalmon", "tomato3",
                       "orange", "darkorange3", "plum1", "purple", "lightgoldenrod", "goldenrod",
                       "darkgreen"), each = 21)
    test <- c("skyblue", "blue4", "darkolivegreen1", "darkolivegreen4","darksalmon", "tomato3",
              "orange", "darkorange3", "plum1", "purple", "lightgoldenrod", "goldenrod",
              "darkgreen")
    
    ## Make plot
    pdf(file.path(currOut_v, currName_v))
    print(cloud(x = y~x+z, 
                data = d, 
                groups = z,
                panel.3d.cloud = panel.3dbars,
                col.facet = test,
                screen = list(z = 75, x = -81, y = 5),
                alpha.facet = .5,
                xbase=0.4, ybase=0.4,
                scales=list(y = list(rot = c(45,45), cex = 0.8),
                            x = list(cex = 0.8),
                            z = list(cex = 1),
                            arrows=FALSE, col=1, distance = c(1.1,.8,.8)),
                xlab = "J Segments",
                ylab = "V Segments",
                zlab = list(label = currLab_v, rot = 90),
                # par.settings = list(axis.line = list(col = "transparent"),
                #                     box.3d = list(col = c(1,1,1,NA,1,NA,1,1,1))),
                par.settings = list(axis.line = list(col = "transparent")),
                main = paste0("Sample: ", currSample_v)))
    dev.off()
    
  } # for i
} # for j


