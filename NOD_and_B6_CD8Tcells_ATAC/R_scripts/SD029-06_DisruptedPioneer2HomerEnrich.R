# goal: pull peaks for motif enrichment from 2 sample LFC testing

library(dplyr)
library(stringr)

motif_db = "HOCOMOCO"
motif_family = read.delim("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_tf-info.txt",
                          header = T)

# ets_terms = c("Ets-like", "Elf-1-like", "Elk-like", "Etv6-like")
# sp1_klf_terms = c("Sp1-like", "Kruppel-like")
# runx_terms = c("Runx")
# bzip_terms = c("Fos factors", "Jun factors", "B-ATF", "ATF-[2-3]-like factors", "NF-E2-like factors")
search_terms = list(ets = "ETS", 
                    #sp1_klf = "SP1/KLF", # no sig SP1/KLF - there was a sign issue with a >/< statement
                    runx = "Runt", bzip = "bZIP")

mean_diff = -0.5 # do we want to impose an LFC for this?
padj = 0.05 #0.1

date.lfc = "2024-09-04"
date.tfdf = "2024-09-03"

phenotype = "LCMV"

consensus_peaks = rtracklayer::import("/Scottbrowne/seq/tmp/devoes/SD029/All_Samples.fwp.filter.non_overlapping.renamed.bed")
# consensus_peaks = rtracklayer::import("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.renamed.bed")
sig.peaks = list.files(path = "/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5",
                       pattern = "all.bed",
                       full.names = T) %>% 
  grep(pattern = "NOD-B6", value = T) %>%
  lapply(., FUN = rtracklayer::import) %>% 
  purrr::reduce(., .f = c) %>% unique()

down.peaks = list.files(path = "/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5",
                       pattern = "down.bed",
                       full.names = T) %>% 
  setNames(nm = basename(.)) %>%
  grep(pattern = "NOD-B6", value = T) %>%
  lapply(FUN = read.delim, header = F, col.names = c("chr","start","end","peakID"))

rm_peak = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/ConsensusPeaks_intersect_AnyHetSNP.bed",
                     header=F)

lfc.res= readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                              motif_db, 
                              "/150bp/padj_0.1/FC_1.5/",
                              phenotype,
                              "/", 
                              date.lfc,
                              "_NoHetSNP_lfc-res_2sample.rds"))
if(phenotype == "Naive"){
  peaks.naive = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/150bp/padj_0.1/FC_1.5/",
                               "Naive/Lost_CART-AND-LCMV_to-Naive_any-strain_peak_subset.rds"))
  
} else{
  peaks.pheno = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/150bp/padj_0.1/FC_1.5/",
                               phenotype,
                               "/Gained_",
                               phenotype,
                               "-to-Naive_any-strain_peak_subset.rds"))
  
}

tf.df = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                            motif_db, 
                            "/150bp/padj_0.1/FC_1.5/",
                            phenotype,
                            "/",
                            date.tfdf,
                            "_tf-df.rds"))
# remove peaks with any HetSNPs
tf.df = tf.df[!(tf.df$GeneID %in% rm_peak$V4),]


# filter lfc with bool
bool_motif = (abs(lfc.res$diff_mean_LFC) > abs(mean_diff)) & (lfc.res$p.adj < padj)
sig_motifs = data.frame(Model = lfc.res$motif[bool_motif])
sig_motifs = left_join(sig_motifs, motif_family, by = "Model")

sites.long = tf.df[, c('GeneID',grep(colnames(tf.df), pattern = "H11MO", value = T))] %>% 
  tidyr::pivot_longer(cols = !GeneID, 
                      names_to = "Model",
                      values_to = "Specificity") %>% 
  rename(peakID="GeneID")

get_peaks <- function(family, sig_motifs, sites.long, lost.peaks, writeBed = F, phenotype = NA, outdir = NA){
  ## define motifs to search for
  search_motifs = sig_motifs %>% filter(Plot_Motif_Color == family)
  # print(head(search_motifs))
  
  ## filter sites for sig motifs in family
  ## filter sites for peaks with B6 specific sites
  sites.keep = sites.long[(sites.long$Model %in% search_motifs$Model) &  sites.long$Specificity == "B6_Specific", ]
  # print(head(sites.keep))
  # print(table(sites.keep$Specificity))
  # print(table(sites.keep$Model))
  sites.keep = sites.keep$peakID %>% unique()
  sites.keep = data.frame(peakID = sites.keep)
  
  ## filter for those with lost signal in NOD
  peaks2homer = inner_join(sites.keep, lost.peaks, by = "peakID")
  
  # write results to bed
  if(writeBed){
    write.table(peaks2homer[,c(2,3,4,1)], 
                file = paste0(outdir, "/",Sys.Date(),"_",phenotype,"_",family %>% str_replace(pattern = "/", replacement = "-"),"_peaks4homer.bed"),
                col.names = F,
                row.names = F,
                quote = F,
                sep = "\t")
  }
  return(peaks2homer)
  
}

peaks4homer = lapply(search_terms,
                     FUN =  get_peaks, 
                     sig_motifs = sig_motifs,
                     sites.long = sites.long,
                     lost.peaks = down.peaks[[grep(phenotype, names(down.peaks), ignore.case = T)]],
                     writeBed = T,
                     outdir = paste0("/Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/homer/",phenotype),
                     phenotype = phenotype)




