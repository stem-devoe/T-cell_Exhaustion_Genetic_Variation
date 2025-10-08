# test LFC ####
testLFC_2sample <- function(motif_name, tf.df, specificity, min_sites, alt = c("less", "greater","two.sided"), parametric = T){
  strain = tf.df[grepl(specificity, tf.df[,motif_name]), ]
  shared = tf.df[grepl("Shared", tf.df[,motif_name]), ]
  # check for minimum number of sites
  if(nrow(strain) < min_sites){
    return(Inf)
    #return(paste0("Too few sites to test: ",nrow(tf.df)))
  }
  
  if(parametric){
    res = t.test(x = strain$logFC, y = shared$logFC, alternative = alt) # default is UNequal
  }else {
    res = wilcox.test(x = strain$logFC, y = shared$logFC, alternative = alt, exact = T)
  }
  
  return(res$p.value)
}

getMeanLFC <- function(motif_name, tf.df, specificity){
  df = data.frame(motif = tf.df[,motif_name],
                  logFC = tf.df[, "logFC"])
  mLFC = df %>% group_by(motif) %>% summarize(mean = mean(logFC))
  if(!any(grepl(specificity, mLFC$motif))){
    return(NA)
  } else{
    return(mLFC$mean[grepl(specificity, mLFC$motif)])
  }
}


