library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(limma)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

set.seed(901)
phenotype = c("CART","LCMV")
common_bg = T
enrich_fill = "LOG_QVALUE" 
enrich_date = "2024-10-05" # date peaks for homer were generated
if(common_bg){
  enrich_path = "/Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/common_bg"
}else{
  enrich_path = "/Scottbrowne/seq/tmp/devoes/SD029/motif_disruption"
}
motif_disruption = c("bZIP",
                     "ETS",
                     "Runt"#,
                     # "SP1-KLF"
)
motif_info = read.delim("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_tf-info.txt")
family_label = setNames(motif_info$Plot_Motif_Color, motif_info$Model)
keep_motif = motif_info$Model[motif_info$keep_enrich_heatmap > 0] # %>% grep(pattern = paste0(c("COT","E2F","Bcl6","Ikzf1","Xbp1","ITF2","STAT1.*1.A"),collapse = "|"), invert = T, value = T)

peaks = list.files(paste0(enrich_path %>% str_remove(pattern = "/common_bg"),"/homer/"),
                   pattern = ".bed",
                   recursive = T,
                   full.names = T) %>% grep(pattern = "Naive", invert = T, value = T)
names(peaks) = basename(peaks) %>% 
  str_remove(pattern = paste0(enrich_date,"_")) %>% 
  str_remove(pattern = "_peaks4.*")
peaks = lapply(peaks, read.delim, header = F, col.names = c("chr","start","end","peakID"))

# upset plot disrupted sites peaks ####

peakIDs = lapply(peaks, FUN = function(peaks){peaks$peakID})

# make combination matrix
comb.mtx = make_comb_mat(peakIDs, mode = "distinct")

# plot
p.upset <- UpSet(comb.mtx,
                 heatmap_height = unit(2.5, "in"),
                 heatmap_width = unit(5, "in"),
                 pt_size = unit(3, "mm"),
                 lwd = 1,
                 set_order = seq_along(set_size(comb.mtx)), # I want the order I manually had
                 comb_order = order(comb_size(comb.mtx), decreasing = T), # do we order that makes sense by connected dots or by most to least 
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize =7
                                     #col = ,
                                     ),
                 top_annotation = upset_top_annotation(comb.mtx, 
                                                       annotation_name_gp = gpar(fontsize = 0),
                                                       annotation_name_rot = 90,
                                                       add_numbers = T,
                                                       axis = F,
                                                       numbers_gp = gpar(fontsize = 6),
                                                       numbers_rot = 90),
                 right_annotation = upset_right_annotation(comb.mtx,
                                                           annotation_name_gp = gpar(fontsize = 0),
                                                          # gp = gpar(fill = ,
                                                          #           col = ),
                                                           add_numbers = T,
                                                           axis = F,
                                                           numbers_gp = gpar(fontsize = 6),
                                                           bar_width = 0.5)
)

p.upset

pdf(file = paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",Sys.Date(),"_MotifDisruptionEnrichment-peaks_Upset.pdf"),
    width = 5,
    height = 3)
p.upset
dev.off()

# scatter ####


enrich_files = list.files(paste0(enrich_path,"/SEA"),
                          pattern = "sea.tsv",
                          full.names = T,
                          recursive = T)  %>% grep(pattern = "Naive",invert = T, value = T)
names(enrich_files) = enrich_files %>% 
  str_remove(pattern = paste0(".*SEA/")) %>%
  str_remove(pattern = "/sea.tsv") %>% str_replace(pattern = "/", replacement = "_")
enrich = lapply(enrich_files, read.delim)

getEnrich <- function(motif_disruption, enrich, enrich_fill, phenotypes){
  pheno1 = getPhenoEnrich(phenotypes[1], motif_disruption, enrich, enrich_fill)
  pheno2 = getPhenoEnrich(phenotypes[2], motif_disruption, enrich, enrich_fill)
  enrich = inner_join(pheno1, pheno2, by = "Model")
  return(enrich)
}

getPhenoEnrich <- function(phenotype, motif_disruption, enrich, enrich_fill){
  enrich = enrich[[paste0(phenotype,"_",motif_disruption)]] %>% filter(ID != "") # account for non tidy data of SEA output
  enrich = enrich[,c('ID',enrich_fill)] 
  colnames(enrich) = c('Model', paste0(phenotype,"_",enrich_fill))
  return(enrich)
}


enrich_scatter = lapply(motif_disruption, FUN = getEnrich, enrich, enrich_fill, phenotype)
names(enrich_scatter) = motif_disruption

enrich_scatter = lapply(enrich_scatter,FUN = left_join, motif_info[,c('Model','Plot_Motif_Color','keep_enrich_heatmap')], by = "Model")

enrich_scatter[[1]] %>% filter(keep_enrich_heatmap > 0 ) %>%
  #mutate(CART_LOG_QVALUE = replace(CART_LOG_QVALUE,CART_LOG_QVALUE < -10, -10)) %>%
  #mutate(LCMV_LOG_QVALUE = replace(LCMV_LOG_QVALUE,LCMV_LOG_QVALUE < -10, -10)) %>%
  ggplot(aes(x = -CART_LOG_QVALUE, y = -LCMV_LOG_QVALUE)) + geom_point() + theme_classic() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') + 
  coord_cartesian(xlim = c(0, 45), 
                                 ylim = c(0, 45)) # is this appropriate? p/q values are very influenced by n of groups. 

                  