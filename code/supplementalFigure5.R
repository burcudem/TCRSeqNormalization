par(mfrow = c(1,1))


library(GGally)
library(MASS)
library(reshape2)
library(viridis)
library(ggplot2)
library(zoo)
library(plyr)
library(dplyr)
library(data.table)

### Data

data_dt <- fread(file.path(baseDir_v, "./data/temp/tempAdaptive.txt"))
outDir_v <- file.path(baseDir_v, "./output/fig6_supfig5/")

### Comparisons



clonality_lsv <- list(c("aClonality", "aClonalityOurCode"),
                      c("aClonality", "ourClonality"),
                      c("aClonalityOurCode", "ourClonality"))

entropy_lsv <- list(c("aEntropy", "aEntropyOurCode"),
                    c("aEntropy", "ourEntropy"),
                    c("aEntropyOurCode", "ourEntropy"))

hyper_lsv <- list(c("aHyper", "ourHyper"))


clonality_lsv_lab<- list(c("Clonality Index - Adaptive by Adaptive","Clonality Index - Adaptive by OTSP"),
                     c("Clonality Index - Adaptive by Adaptive","Clonality Index - OTSP by OTSP"),
                     c("Clonality Index - Adaptive by OTSP","Clonality Index - OTSP by OTSP"))

entropy_lsv_lab <- list(c("Shannon Diversity Index - Adaptive by Adaptive","Shannon Diversity Index - Adaptive by OTSP"),
                        c("Shannon Diversity Index - Adaptive by Adaptive","Shannon Diversity Index - OTSP by OTSP"),
                        c("Shannon Diversity Index - Adaptive by OTSP","Shannon Diversity Index - OTSP by OTSP"))
hyper_lsv_lab <- list(c("Hyperexpanded Clones - Adaptive Platform", "Hyperexpanded Clones - OTSP"))

compare_lslsv <- list("Clonality" = clonality_lsv,
                      "Entropy" = entropy_lsv,
                      "Hyper" = hyper_lsv)

compare_lslsv_lab <- list("Clonality" = clonality_lsv_lab,
                      "Entropy" = entropy_lsv_lab,
                      "Hyper" = hyper_lsv_lab)

for (i in 1:length(compare_lslsv)) {
  
  ## Get name and list
  currMetric_v <- names(compare_lslsv)[i]
  currList_lsv <- compare_lslsv[[currMetric_v]]
  currList_lsv_lab<-compare_lslsv_lab[[currMetric_v]]
  ## Make directory
  #mkdir <- p /path-to-directory/directory-name#
  

    currDir_v <- file.path(outDir_v,currMetric_v) 
    print(currDir_v)
    dir.create(currDir_v, recursive = TRUE)
    
    for (j in 1:length(currList_lsv)) {
    
    ## Get columns
    currCol1_v <- currList_lsv_lab[[j]][1]
    currCol2_v <- currList_lsv_lab[[j]][2]
    ## Get data
    currVals1_v <- data_dt[[currList_lsv[[j]][1]]]
    currVals2_v <- data_dt[[currList_lsv[[j]][2]]]
    
    ## Get correlation
    currCor <- cor.test(currVals1_v, currVals2_v)
    currlm <- summary(lm(currVals2_v ~ currVals1_v))
   
    ## Get results
    currP_v <- currCor$p.value
    currR_v <- currCor$estimate
    
    ## Get line paramters
    currA_v <- currlm$coefficients[1,1]
    currB_v <- currlm$coefficients[2,1]
    
    ## Make legend
    currLegend_v <- c(paste0("R = ", round(currR_v, digits = 1)))
                      #paste0("p = ", round(currP_v, digits = 1)))
    
    
    ## Make plot info
    title_v <- paste("OTSP vs Adaptive")
    name_v <- file.path(currDir_v, paste0(currCol1_v, "_vs_", currCol2_v, ".pdf"))
    
    ## Make plot
    pdf(file = name_v)
    plot(x = currVals1_v, xlab = currCol1_v,
         y = currVals2_v, ylab = currCol2_v,
         main = title_v, cex.lab=1.3, cex.main=1.3)
    legend("bottomright", legend = currLegend_v, bty = 'n', text.font=2, text.col = "Red" , c(as.expression(bquote(bold("Bold")))))
    abline(a = currA_v, b = currB_v)
    dev.off()
  } # for j
} # for i



