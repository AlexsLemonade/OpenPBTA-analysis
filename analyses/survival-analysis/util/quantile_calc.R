# Functions for calculating quantiles of variables used for survival analysis

# using 0.5 as a cutoff
calc_groups <- function(meta_indep, tp53_quantiles, telomerase_quantiles){
  meta_indep$tp53_0.5 <- ifelse(meta_indep$tp53_score > 0.5, "tp53 high",
                                ifelse(meta_indep$tp53_score <= 0.5, "tp53 low",NA))
  
  # using 0.5 as a cutoff
  meta_indep$tel_0.5 <- ifelse(meta_indep$telomerase_score > 0.5, "telomerase high",
                               ifelse(meta_indep$telomerase_score <= 0.5, "telomerase low",NA))
  
  # tp53 strict
  meta_indep$tp53_strict <- ifelse(meta_indep$tp53_score >= tp53_quantiles[[4]], "tp53_high",
                                   ifelse(meta_indep$tp53_score > tp53_quantiles[[2]] & meta_indep$tp53_score < tp53_quantiles[[4]], "tp53_mid",
                                          ifelse(meta_indep$tp53_score <= tp53_quantiles[[2]], "tp53_low", NA)))
  table(meta_indep$tp53_strict)
  
  # tel strict
  meta_indep$tel_strict <- ifelse(meta_indep$telomerase_score >= telomerase_quantiles[[4]], "telomerase_high",
                                  ifelse(meta_indep$telomerase_score >  telomerase_quantiles[[2]] & meta_indep$telomerase_score < telomerase_quantiles[[4]], "telomerase_mid",
                                         ifelse(meta_indep$telomerase_score <= telomerase_quantiles[[2]], "telomerase_low", NA)))
  table(meta_indep$tel_strict)
  
  
  # using 0.5 as a cutoff
  meta_indep$pheno_0.5 <- ifelse(meta_indep$tp53_score > 0.5 & meta_indep$telomerase_score > 0.5, "tp53_high telomerase_high",
                                 ifelse(meta_indep$tp53_score > 0.5 &  meta_indep$telomerase_score <= 0.5, "tp53_high telomerase_low",
                                        ifelse(meta_indep$tp53_score <= 0.5 &  meta_indep$telomerase_score > 0.5, "tp53_low telomerase_high",
                                               ifelse(meta_indep$tp53_score <= 0.5 &  meta_indep$telomerase_score <= 0.5, "tp53_low telomerase_low",
                                                      NA))))
  table(meta_indep$pheno_0.5)
  
  # using strict cutoffs
  meta_indep$pheno_strict <- ifelse(meta_indep$tp53_score >= tp53_quantiles[[4]] & meta_indep$telomerase_score >= telomerase_quantiles[[4]], "tp53_high telomerase_high",
                                    ifelse(meta_indep$tp53_score >= tp53_quantiles[[4]] &  meta_indep$telomerase_score <= telomerase_quantiles[[2]], "tp53_high telomerase_low",
                                           ifelse(meta_indep$tp53_score <= tp53_quantiles[[2]] &  meta_indep$telomerase_score >= telomerase_quantiles[[4]], "tp53_low telomerase_high",
                                                  ifelse(meta_indep$tp53_score <= tp53_quantiles[[2]] &  meta_indep$telomerase_score <= telomerase_quantiles[[2]], "tp53_low telomerase_low",
                                                         ifelse(meta_indep$tp53_score > tp53_quantiles[[2]] & meta_indep$tp53_score < tp53_quantiles[[4]] &  meta_indep$telomerase_score <= telomerase_quantiles[[2]], "tp53_mid telomerase_low",
                                                                ifelse(meta_indep$tp53_score <= tp53_quantiles[[2]] &  meta_indep$telomerase_score > telomerase_quantiles[[2]] & meta_indep$telomerase_score < telomerase_quantiles[[4]], "tp53_low telomerase_mid",
                                                                       "tp53 tel mid"))))))
  
  table(meta_indep$pheno_strict)
  
  ### without mid
  meta_indep$pheno_extremes <- ifelse(meta_indep$tp53_score >= tp53_quantiles[[4]] & meta_indep$telomerase_score >= telomerase_quantiles[[4]], "tp53_high telomerase_high",
                                      ifelse(meta_indep$tp53_score >= tp53_quantiles[[4]] &meta_indep$telomerase_score <= telomerase_quantiles[[2]], "tp53_high telomerase_low",
                                             ifelse(meta_indep$tp53_score <= tp53_quantiles[[3]] & meta_indep$telomerase_score >= telomerase_quantiles[[4]], "tp53_low telomerase_high",
                                                    ifelse(meta_indep$tp53_score <= tp53_quantiles[[2]] & meta_indep$telomerase_score <= telomerase_quantiles[[2]], "tp53_low telomerase_low", NA))))
  
  ### get the results out 
  return(meta_indep)
}
