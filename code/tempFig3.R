### Holding onto some old checks just in case

### Check counts
myNew <- data_lsdt$LIB170111LC_new_dil
my5 <- data_lsdt$LIB170111LC_dil5

myNew$VJ <- paste(gsub("-", "", myNew$V), myNew$J, sep = "_")
my5$VJ <- paste(gsub("-", "", my5$V), my5$J, sep = "_")

checkNewCounts <- merge(myNew, newDil, by = "VJ")
check5Counts <- merge(my5, dil5, by = "VJ")

for (i in 125:134) {
  mine_v <- paste0("S", i)
  old_v <- as.character(i)
  check_v <- check5Counts[[mine_v]] == check5Counts[[old_v]]
  bad_v <- which(!check_v)
  if (length(bad_v) > 0) {
    print(i)
    print(bad_v)
  }
}


### Check estimates
checkOld <- merge(combined[,c("VJ", "mu.nb_old.Estimate", "theta_old")],
                  disp_df[,c("VJ", "mu.nb_unknown.Estimate", "theta_unknown")], by = "VJ")
which(checkOld$mu.nb_old.Estimate != checkOld$mu.nb_unknown.Estimate)
which(checkOld$theta_old != checkOld$theta_unknown)
checkDil5 <- merge(combined[,c("VJ", "mu.nb_dil5.Estimate", "theta_dil5")],
                   disp_df[,c("VJ", "mu.nb_dil5.Estimate", "theta_dil5")], by = "VJ")
which(checkDil5$mu.nb_dil5.Estimate.x != checkDil5$mu.nb_dil5.Estimate.y)
which(checkDil5$theta_dil5.x != checkDil5$theta_dil5.y)
checkNew <- merge(combined[,c("VJ", "mu.nb_newdil.Estimate", "theta_newdil")],
                  disp_df[,c("VJ", "mu.nb_new_dil.Estimate", "theta_new_dil")], by = "VJ")
which(checkNew$mu.nb_newdil.Estimate != checkNew$mu.nb_new_dil.Estimate)
which(checkNew$theta_newdil != checkNew$theta_new_dil)