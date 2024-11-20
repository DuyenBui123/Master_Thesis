DRMAT_conmax <- function(bpTfile, bpSfile) {
  bp_S <- read.delim(bpSfile,header = F)
  
  bp_T <- read.delim(bpTfile,header = F)
  
  bp_T_new <- bp_T %>% pivot_longer(!V1, names_to = "col", values_to = "breakpoint")
  bp_T_new <- na.omit(bp_T_new)
  bp_T_new <- subset(bp_T_new, select = -c(col))
  colnames(bp_T_new) <- c("sample_id", "Breakpoint")
  bp_T_new <- bp_T_new[!bp_T_new$Breakpoint == -999,]
  
  bp_S_new <- bp_S %>% pivot_longer(!V1, names_to = "col", values_to = "breakpoint")
  bp_S_new <- subset(na.omit(bp_S_new), select = -c(col))
  colnames(bp_S_new) <- c("sample_id", "Breakpoint")
  bp_S_new <- bp_S_new[!bp_S_new$Breakpoint == -999,]
  
  bp_IC_ridgidnr <- rbind(bp_T_new, bp_S_new)
  
  if (length(unique(cal_comb_ref_rm$sample_id)) == nrow(cal_ndvi_nonNan)) {
    DRMAT_SR_cal <- validateAlgorithmTotal(ref_df = cal_comb_ref_rm, 
                                           algo_df = bp_IC_ridgidnr, cl= mycores)
    DRMAT_IC_ridgidnr_stats <- myFPStats(DRMAT_SR_cal, NewAccuracy = TRUE)
    
  }
  return(DRMAT_IC_ridgidnr_stats)
}