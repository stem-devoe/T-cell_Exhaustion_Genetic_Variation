# categorize peaks by signal change ####

categorizePeakChange <- function(decideTests_row){
  # print(decideTests_row)
  # c(a,b,c,...) == c(x,y,z,..) will evaluate each one (a==x, b==y, c==z)
  case_when(
    all(decideTests_row == c(1,1)) ~ "GG",
    all(decideTests_row == c(1,0)) ~ "GN",
    all(decideTests_row == c(0,1)) ~ "NG",
    all(decideTests_row == c(-1,-1)) ~ "LL",
    all(decideTests_row == c(-1,1)) ~ "LG",
    all(decideTests_row == c(1,-1)) ~ "GL",
    all(decideTests_row == c(-1,0)) ~ "LN",
    all(decideTests_row == c(0,-1)) ~ "NL",
    sum(abs(decideTests_row)) == 0  ~ "NN"
  )
}

# categorize peaks by sites ####
readMotifSites <- function(motif_path, motif_db, motif_name){
  motif_gr = read.delim(paste0(motif_path,
                               "/ConsensusPeaks_",motif_db,"_",motif_name,"_categorized.bed"), # reusability issue!!!
                        header = F,
                        col.names = c("chr","start","end","originalCoordName",
                                      "strand","motif","Specificity")) 
  motif_gr = makeGRangesFromDataFrame(df = motif_gr[,c(1:4,7)], 
                                      ignore.strand = T, 
                                      starts.in.df.are.0based = T,
                                      keep.extra.columns = T)
  
  return(motif_gr)
}

intersectPeak2Motif <- function(motif_gr, peakset_gr){
  # split motif sites by specificity
  motif_gr = split(motif_gr, f = as.factor(motif_gr$Specificity))
  
  # intersect with peaks
  overlaps = lapply(motif_gr, 
                    FUN = findOverlaps,
                    subject = peakset_gr,
                    type = "any",
                    select = "all") %>% lapply(X = .,
                                               FUN = function(x, peak_gr){peak_gr$name[subjectHits(x)]},
                                               peak_gr = peakset_gr)
  
  return(overlaps)
  
}

matchPeak2Site <- function(site_specificity, overlaps_res, peak_names){
  site = rep(0, length(peak_names))
  #print(head(site))
  names(site) = peak_names
  #print(head(site))
  #print(head(overlaps.res[[site.specificity]]))
  site[names(site) %in% unique(overlaps_res[[site_specificity]])] = 1
  # print(length(site))
  return(site)
}

formatPeakSites <- function(overlaps, peak_sites, peakset_gr){
  names(peak_sites) = names(overlaps)
  peak_sites = data.frame(peak_sites)
  # fill in missing specificities
  if(ncol(peak_sites) < 4){
    specs = colnames(peak_sites)
    missing.specs = setdiff(c("B6_Specific","Common","Common_Polymorphic","NOD_Specific"),
                            specs)
    for(i in seq_along(missing.specs)){
      peak_sites[,missing.specs[i]] = 0
    }
  }
  # order
  peak_sites = peak_sites[,c("B6_Specific","Common","Common_Polymorphic","NOD_Specific")]
  
  rownames(peak_sites) = peakset_gr$name
  return(peak_sites)
}

assignSiteCategory <- function(sites_row){
  # sites_row must be in order of B6_Specific, Common, Common_Polymorhpic, NOD_Specific
  # c(a,b,c,...) == c(x,y,z,..) will evaluate each one (a==x, b==y, c==z)
  case_when(
    # B6-specific and NOD-specific site within same peak
    sites_row[1] & sites_row[4] ~ "Mixed_Specificity",
    # B6-specific site
    sites_row[1] & !sites_row[4] ~ "B6_Specific",
    # NOD-specific site
    !sites_row[1] & sites_row[4] ~ "NOD_Specific",
    # Common but polymorphic
    #sites_row[3] == 1 ~ "Common_Polymorphic", 
    sites_row[3] == 1 ~ "Shared", 
    # No Site
    sum(abs(sites_row)) == 0  ~ 'No_Site',
    # Common
    #sites_row[2] == 1 ~ "Common",
    sites_row[2] == 1 ~ "Shared", #"Common"
    TRUE ~ "Failed"
  )
}

getSingleMotifCategories <- function(motif_name, motif_db, motif_path, peakset_gr){
  
  motif_gr = readMotifSites(motif_path, motif_db, motif_name)
  
  overlaps = intersectPeak2Motif(motif_gr, peakset_gr)
  
  peak_sites = lapply(names(overlaps),
                      FUN = matchPeak2Site,
                      overlaps_res = overlaps, 
                      peak_names = peakset_gr$name)
  
  peak_sites = formatPeakSites(overlaps, peak_sites, peakset_gr)
  
  # return(peak_sites)
  
  peak_specificity = apply(peak_sites, 1, FUN = assignSiteCategory) 
  
  peak_spec.df = data.frame(motif_spec = peak_specificity,
                            GeneID = rownames(peak_sites))
  colnames(peak_spec.df)[1] = motif_name
  return(peak_spec.df)
}

getMotifCategories <- function(motif_names, motif_db, motif_path, peakset_gr, ncores = NULL){
  
  if(is.null(ncores)){
    motif_mtx.list = lapply(X = motif_names,
                            FUN = getSingleMotifCategories,
                            motif_db,
                            motif_path,
                            peakset_gr
    )
  } else {
    motif_mtx.list = parallel::mclapply(X = motif_names,
                                        FUN = getSingleMotifCategories,
                                        motif_db,
                                        motif_path,
                                        peakset_gr,
                                        mc.cores = ncores
    )
  }
  
  motif_mtx = purrr::reduce(motif_mtx.list, left_join, by = "GeneID")
  rownames(motif_mtx) = motif_mtx$GeneID
  return(motif_mtx)
}


# test LFC ####
testLFC <- function(motif_name, tf.df, specificity, min_sites, alt = c("less", "greater","two.sided"), parametric = T){
  tf.df = tf.df[grepl(specificity, tf.df[,motif_name]), ]
  # check for minimum number of sites
  if(nrow(tf.df) < min_sites){
    return(Inf)
    #return(paste0("Too few sites to test: ",nrow(tf.df)))
  }
  
  if(parametric){
    res = t.test(tf.df$logFC, alternative = alt)
  }else {
    res = wilcox.test(tf.df$logFC, alternative = alt, exact = T)
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


# Plotting ####

plotTFA <- function(motif_name, tf.df, tfa.res){
  
  tf.df = tf.df[!grepl("Mixed|No_Site", tf.df[,motif_name]), c('logFC', motif_name)]
  colnames(tf.df) = c("logFC","motif_specificity")
  ggplot(tf.df, aes(x = motif_specificity, fill = motif_specificity, y = logFC)) + 
    geom_boxplot() + 
    geom_jitter(position = position_dodge(width = 0.75)) +
    ggtitle(motif_name,
            subtitle = paste0("B6_Specific padj = ", tfa.res$p.adj[tfa.res$motif %in% motif_name])) +
    theme_classic()
}
