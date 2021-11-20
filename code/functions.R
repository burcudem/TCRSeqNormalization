################################################################################
################################## FUNCTIONS ###################################
################################################################################
### 1. mergeDTs      - merge list of data.tables                             ###
### 1. convertDFT    - convert data.frame/matrix to data.table or vice versa ###
### 2. big_label     - theme for ggplot                                      ###
### 4. my_theme      - theme for ggplot                                      ###
### 5. calcDisp      - calculate NB parameters                               ###
################################################################################
################################################################################
################################################################################

###
### MERGEDTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

mergeDTs <- function(data_lsdt, mergeCol_v, keepCol_v = NULL, ...) {
  #' Merge many data.tables together
  #' @description Take many data.tables and merge on and ID column, 
  #' extracting a single column from each data.table as the column of interest
  #' @param data_lsdt list of data.tables to merge
  #' @param mergeCol_v which column from all of the data.tables to use to merge
  #' @param keepCol_v which column from all of the data.tables to use as the column of interest. If NULL, use all columns
  #' @param ... extra parameters passed to merge
  #' @return data.table with ncol == length(data_lsdt) + 1. Column names are names of list, or defaults to V1, V2,...
  #' @export
  
  ## Grab extra arguments
  extraParams_lsv <- list(...)
  
  ## Handle extra arguments
  if (!is.null(extraParams_lsv$all)){
    all_v <- extraParams_lsv$all
  } else {
    all_v <- T
  } # fi
  
  if (!is.null(extraParams_lsv$sort)){
    sort_v <- extraParams_lsv$sort
  } else {
    sort_v <- F
  } # fi
  
  ## If keepCol_v is NULL, grab all other columns
  if (is.null(keepCol_v)){
    keepCol_v <- colnames(data_lsdt[[1]])[-which(colnames(data_lsdt[[1]]) %in% mergeCol_v)]
  } # fi
  
  ## Create initial table by extracting the 2 columns of interest from the rest
  merge_dt <- data_lsdt[[1]][,mget(c(mergeCol_v, keepCol_v))]
  
  ## Create initial column names (first check if list has names and add if not)
  if (is.null(names(data_lsdt))) {
    names_v <- paste("V", 1:length(data_lsdt))
    names(data_lsdt) <- names_v
  } # fi
  
  if (length(keepCol_v) > 1){
    colNames_v <- c(mergeCol_v, paste(names(data_lsdt)[1], keepCol_v, sep = "_"))
  } else {
    colNames_v <- c(mergeCol_v, names(data_lsdt)[1])
  } # fi
  
  for (i in 2:length(data_lsdt)) {
    
    ## This is new (2018-10-10) - need to make new keepCol_v if the data.tables don't have same columns
    if (!keepCol_v %in% colnames(data_lsdt[[i]])) {
      keepCol_v <- colnames(data_lsdt[[i]])[-which(colnames(data_lsdt[[i]]) %in% mergeCol_v)]
    } # fi
    
    ## Merge
    merge_dt <- merge(merge_dt,
                      data_lsdt[[i]][,mget(c(mergeCol_v, keepCol_v))],
                      by = mergeCol_v,
                      all = all_v, sort = sort_v)
    ## Update column names
    if (length(keepCol_v) > 1){
      colNames_v <- c(colNames_v, paste(names(data_lsdt)[i], keepCol_v, sep = "_"))
    } else {
      colNames_v <- c(colNames_v, names(data_lsdt)[i])
    } # fi
    
    ## Rename columns
    colnames(merge_dt) <- colNames_v
  } # for i
  return(merge_dt)
} # mergeDTs

###
### CONVERTDFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

convertDFT <- function(data_dft, col_v = NA, newName_v = "V1", rmCol_v = T) {
  #' Convert between data.table and data.frame
  #' @description 
  #' Change data.tables into data.frames with specified row.names or
  #' data.frames/matrices into data.tables, copying over row.names as the first column.
  #' If data.frame/matrix doesn't have row.names, then no columns will be added to data.table
  #' @param data_dft data in either data.table or data.frame format (can also be a matrix)
  #' @param col_v character or numeric vector. if converting from dt to df, column name or index of which column to use as row.names.
  #' NA (default) will use 1st column; NULL will not add rownames
  #' @param newName_v character vector. if converting from df/mat to dt, what to name new column. (default is "V1")
  #' if newName_v is already a column name, will paste "_2" to end of newName_v.
  #' @param rmCol_v boolean value indicating whether to remove the column used to make the rownames from the output table (T) or to leave it (F)
  #' @return either a data.table or data.frame (opposite class of input)
  #' @examples 
  #' # Data
  #' my_df <- data.frame("A" = 1:10, "B" = LETTERS[1:10], "C" = letters[11:20])
  #' my_df2 <- my_df; rownames(my_df2) <- paste0("Row", 1:10)
  #' my_mat <- as.matrix(my_df); my_mat2 <- as.matrix(my_df2)
  #' my_dt <- data.table("AA" = 10:1, "BB" = LETTERS[5:14], "CC" = letters[20:11])
  #' convertDFT(data_dft = my_df)
  #' convertDFT(data_dft = my_df2)
  #' convertDFT(data_dft = my_df2, newName_v = "Test")
  #' convertDFT(data_dft = my_mat2, newName_v = "MatTest")
  #' convertDFT(data_dft = my_dt)
  #' convertDFT(data_dft = my_dt, col_v = NULL)
  #' convertDFT(data_dft = my_dt, col_v = "BB")
  #' convertDFT(data_dft = my_dt, col_v = 3, rmCol_v = F)
  #' @export
  
  ## Row names function
  addRowNames <- function(data_dft, out_dft, newName_v) {
    newName_v <- ifelse(newName_v %in% colnames(out_dft), paste0(newName_v, "_2"), newName_v)
    out_dft[[newName_v]] <- rownames(data_dft)
    out_dft <- out_dft[, c(ncol(out_dft), 1:(ncol(out_dft)-1)), with = F]
  } # addRowNames
  
  ## Get class
  class_v <- class(data_dft)
  
  ## Convert data.table to data.frame
  if ("data.table" %in% class_v){
    
    ## Convert
    out_dft <- as.data.frame(data_dft)
    
    ## Get column for row names
    col_v <- ifelse(is.na(col_v), colnames(data_dft)[1],
                    ifelse(is.null(col_v), NULL,
                           ifelse(is.numeric(col_v), colnames(data_dft)[col_v], col_v)))
    
    ## Add row names and handle column that provided names
    if (length(col_v) > 0) {
      
      rownames(out_dft) <- data_dft[[col_v]]
      
      ## Remove column that provided rownames
      if (rmCol_v) {
        whichCol_v <- which(colnames(data_dft) == col_v)
        out_dft <- out_dft[,-whichCol_v, drop = F]
      } # fi
      
    } # fi
    
    ## Convert data.frame to data.table
  } else if (class_v == "data.frame"){
    
    ## Convert
    out_dft <- as.data.table(data_dft)
    
    ## Handle row names
    if (!identical(rownames(data_dft), as.character(1:nrow(data_dft)))) {
      out_dft <- addRowNames(data_dft, out_dft, newName_v)
    } # fi
    
    ## Convert matrix to data.table
  } else if (class_v == "matrix") {
    
    ## Convert
    out_dft <- as.data.table(data_dft)
    
    ## Handle row names
    if (!is.null(rownames(data_dft))) {
      out_dft <- addRowNames(data_dft, out_dft, newName_v)
    } # fi
    
  } else {
    stop("Neither 'data.table', 'data.frame', nor 'matrix' were in the class of data_dft. Please check your input data.")
  } # fi
  
  ## Return
  return(out_dft)
  
} # convertDFT


###
### BIGLABEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

big_label <- function() {
  #' Big Label Theme
  #' @description Same as my_theme(), but even larger text and also y-axis labels are angled 45
  #' @export
  
  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(size = 16),
          axis.text.y = element_text(angle = 45),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
} # big_label

###
### MYTHEME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

my_theme <- function() {
  #' Base Custom Theme
  #' @description Minor changes to theme_classic() - slightly larger text and centered title
  #' @export
  
  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
}

###
### CALCDISP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

calcDisp <- function(spike, names) {
  
  ## Calculations
  m.nb<-glm.nb (spike ~ 1)
  theta<-(summary(m.nb)$theta)
  mu.nb<- coef(summary(m.nb, dispersion=1)) 
  
  ## Turn to data.frames
  mu_df <- data.frame(mu.nb)
  colnames(mu_df) <- paste(names[1], colnames(mu_df), sep = '.')
  theta_df <- data.frame(theta)
  colnames(theta_df) <- names[2]
  
  ## Make results
  results <- cbind(mu_df, theta_df)
  
  ## Output
  return(results)
} # calcDisp