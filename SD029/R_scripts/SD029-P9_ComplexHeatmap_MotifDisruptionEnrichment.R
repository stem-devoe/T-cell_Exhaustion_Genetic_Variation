library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(limma)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

set.seed(901)
phenotype = "CART"
center_signal = T
common_bg = T
enrich_fill = "LOG_QVALUE" # "ENR_RATIO"
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

v = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/v_150bp.rds")
norm_counts = v$E %>% as.data.frame()
colnames(norm_counts) = colnames(norm_counts) %>% str_remove(pattern = "_REP1.*")

meta = read.delim("/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv")
groups_key = setNames(meta$Group, meta$SampleID)

peaks = list.files(paste0(enrich_path %>% str_remove(pattern = "/common_bg"),"/homer/",phenotype),
                   pattern = ".bed",
                   full.names = T)
names(peaks) = basename(peaks) %>% 
  str_remove(pattern = paste0(enrich_date,"_",phenotype,"_")) %>% 
  str_remove(pattern = "_peaks4.*")
peaks = lapply(peaks, read.delim, header = F, col.names = c("chr","start","end","peakID"))

# signal heatmap ####

# (may have peaks represented >1 time)

getSignal <- function(motif_disruption, peaks, signal, groups_key, filter_groups = F, keep_groups = NULL){
  peaks_keep = peaks[[motif_disruption]]$peakID
  signal = signal[peaks_keep,] %>% mutate(peakID = rownames(.))
  
  if(filter_groups){
    signal = signal[,c(names(groups_key[grepl(paste0(keep_groups,collapse="|"),groups_key)]),"peakID")]
  }
  
  
  df = tidyr::pivot_longer(signal,
                           cols = !peakID,
                           names_to = "SampleID",
                           values_to = "Signal") %>%
    mutate(Group = groups_key[SampleID], motif_disruption = motif_disruption)
  
  return(df)
  
}

signal_long = lapply(motif_disruption, getSignal, peaks, norm_counts, groups_key, 
                     filter_groups = T, 
                     keep_groups = c(phenotype,"Naive")) %>%
  purrr::reduce(rbind)
signal_long.mean = signal_long %>% group_by(Group, peakID, motif_disruption) %>% summarise(mean = mean(Signal)) %>% mutate(ploty = paste0(motif_disruption,"_",peakID))
signal_wide.mean = signal_long.mean %>% 
  tidyr::pivot_wider(names_from = Group, values_from = mean) %>% 
  relocate(B6_Naive, .before = NOD_Naive)

signal_mtx = as.matrix(signal_wide.mean[,grepl(paste0(phenotype,"|Naive"),colnames(signal_wide.mean))])
rownames(signal_mtx) = signal_wide.mean$ploty

if(center_signal){
  signal_mtx = signal_mtx %>% t() %>% scale( center = T, scale = F) %>% t()
}

hm.signal = Heatmap(signal_mtx,
                    # remove clustering and dendrograms
                    cluster_rows = F,
                    show_row_names = F,
                    cluster_columns = F,
                    # set order of peaks to group by disruption
                    row_order = order(rownames(signal_mtx)),
                    # 
                    row_split = signal_wide.mean$motif_disruption,
                    gap = unit(2,"mm"),
                    border = "black",
                    height = unit(2, "in"), # 4
                    col = circlize::colorRamp2(breaks = c(min(signal_mtx),0,max(signal_mtx)), c("blue","white","red"), reverse = FALSE),
                    width = unit(0.15 * ncol(signal_mtx), "in") # 0.25
)
hm.signal

# logFC plot ####

tfit.strain = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/150bp/tfit.rds")
tt.strain = topTreat(tfit.strain, coef = paste0('NOD-B6_',phenotype), n = Inf) %>% 
  rename(peakID = GeneID)

lfc = signal_long.mean %>% 
  left_join(y = tt.strain[,c('peakID','logFC')], by = "peakID") %>% 
  ungroup() %>% 
  select(!Group & !mean) %>% 
  unique %>% 
  mutate(plotx = "NOD vs B6") 

hm.lfc = Heatmap(as.matrix(lfc[,'logFC']),
                 # remove clustering and dendrograms
                 cluster_rows = F,
                 show_row_names = F,
                 cluster_columns = F,
                 # set order of peaks to group by disruption
                 row_order = order(lfc$ploty), # same labels as signal so should order the same
                 # 
                 row_split = lfc$motif_disruption,
                 gap = unit(2,"mm"),
                 border = "black",      
                 col = circlize::colorRamp2(breaks = c(min(lfc$logFC),0,0.1), c("salmon4","white","darkcyan"), reverse = FALSE),
                 height = unit(2, "in"), # 4
                 width = unit(0.15, "in") # 0.25
)
hm.lfc

s# enrichment heatmap ####

enrich_files = list.files(paste0(enrich_path,"/SEA/",phenotype),
                          pattern = "sea.tsv",
                          full.names = T,
                          recursive = T) 
names(enrich_files) = enrich_files %>% 
  str_remove(pattern = paste0(".*",phenotype,"/")) %>%
  str_remove(pattern = "/sea.tsv")
enrich = lapply(enrich_files, read.delim)

getEnrich <- function(motif_disruption, enrich, enrich_fill){
  enrich = enrich[[motif_disruption]] %>% filter(ID != "") # account for non tidy data of SEA output
  enrich = enrich[,c('ID',enrich_fill)]
  enrich$motif_disruption = motif_disruption
  return(enrich)
}

enrich_long = lapply(motif_disruption, getEnrich, enrich, enrich_fill)  %>% purrr::reduce(rbind)
enrich_long = enrich_long[enrich_long$ID %in% keep_motif, ]  %>% arrange(ID)

if(enrich_fill == "LOG_QVALUE"){
  enrich_wide = enrich_long %>% 
    mutate(LOG_QVALUE = replace(LOG_QVALUE, LOG_QVALUE < -10, -10)) %>% 
    tidyr::pivot_wider(names_from = ID, values_from = LOG_QVALUE) 
} else if(enrich_fill == "ENR_RATIO"){
  enrich_wide = enrich_long %>% 
    mutate(ENR_RATIO = replace(ENR_RATIO, ENR_RATIO > 10, 10)) %>% 
    tidyr::pivot_wider(names_from = ID, values_from = ENR_RATIO) 
  
}

# duplicate rows to match num peaks per disrupted motif to get uniform split sizes
num_rows = table(signal_wide.mean$motif_disruption)
for(i in seq_along(num_rows)){
  j = 1
  feature = names(num_rows)[i]
  while(j < num_rows[i]){
    enrich_wide = rbind(enrich_wide, enrich_wide[grep(feature,enrich_wide$motif_disruption)[1],] )
    j = j + 1
  }
}
table(enrich_wide$motif_disruption)

enrich_mtx = enrich_wide[,grep("H11MO",colnames(enrich_wide))] %>% as.matrix()
rownames(enrich_mtx) = enrich_wide$motif_disruption

if(enrich_fill == "LOG_QVALUE"){
  enrich_col = circlize::colorRamp2(breaks = c(0,-1*min(enrich_mtx)), c("yellow","darkgreen"), reverse = FALSE)
  enrich_mtx = -1 * enrich_mtx
}else if(enrich_fill == "ENR_RATIO"){
  enrich_col = circlize::colorRamp2(breaks = c(0,max(enrich_mtx)), c("yellow","darkgreen"), reverse = FALSE)
}


hm.enrich = Heatmap(enrich_mtx,
        # remove clustering and dendrograms
        cluster_rows = F,
        show_row_names = F,
        cluster_columns = F,
        # set order of peaks to group by disruption
        row_order = order(enrich_wide$motif_disruption), # same labels as signal so should order the same
        # 
        row_split = enrich_wide$motif_disruption,
        gap = unit(2,"mm"),
        border = "black",    
        
        column_order = order(family_label[colnames(enrich_mtx)]),
        column_split = family_label[colnames(enrich_mtx)],
        column_title_gp = gpar(fontsize = 8),
        column_title_rot = 90,
        #column_names_rot = 45,
        column_names_gp = gpar(fontsize = 6), # 8
        
        col = enrich_col,
        height = unit(2,"in"), # 4
        width = unit(5, #0.075*length(keep_motif),
                     "in")  # 0.1
)
hm.enrich

# enrichment dotplot ####

getEnrichDotPlot <- function(motif_disruption, enrich){
  enrich = enrich[[motif_disruption]] %>% filter(ID != "") # account for non tidy data of SEA output
  enrich = enrich[,c('ID','QVALUE','LOG_QVALUE','ENR_RATIO')]
  enrich$motif_disruption = motif_disruption
  return(enrich)
}

enrich_dp = lapply(motif_disruption, getEnrichDotPlot, enrich)  %>% purrr::reduce(rbind)
enrich_dp = enrich_dp[enrich_dp$ID %in% keep_motif, ]  %>%
  rename(Model = "ID") %>%
  left_join(y = motif_info[,c('Model','Plot_Motif_Color')], by = "Model") %>% arrange(Plot_Motif_Color,Model) %>%
  mutate(motif_disruption = factor(motif_disruption, levels = c("Runt","ETS","bZIP")),
         Model = factor(Model, levels = unique(Model))) 

dp.enrich = enrich_dp %>%
  mutate(LOG_QVALUE = replace(LOG_QVALUE, LOG_QVALUE < -5, -5)) %>% 
  mutate(ENR_RATIO = replace(ENR_RATIO, ENR_RATIO > 10, 10)) %>%
  ggplot(aes(x = Model, y = motif_disruption, size = ENR_RATIO, color = -LOG_QVALUE)) +
  geom_point() + 
  scale_color_viridis_c() +
  theme_classic() + my.theme +
  theme(axis.text.x = element_text(angle = 90)) # facet on motif family?
dp.enrich

dp.enrich.egg = egg::set_panel_size(dp.enrich, 
                                    height = unit(4,"in"),
                                    width = unit(0.1*length(keep_motif),"in"))

# save plots ####

# ComplexHeatmap Does not Work with egg
if(center_signal){
  signal_file = paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",
                       Sys.Date(),
                       "_",
                       phenotype,
                       "_hm-centered-signal_motif-disruption.pdf")
}else{
  signal_file = paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",
                    Sys.Date(),
                    "_",
                    phenotype,
                    "_hm-signal_motif-disruption.pdf")
}

pdf(file = signal_file,
    width = 8.5,
    height = 11)
hm.signal
dev.off()


pdf(file = paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",
                  Sys.Date(),
                  "_",
                  phenotype,
                  "_hm-lfc_motif-disruption.pdf"),
    width = 8.5,
    height = 11)
hm.lfc
dev.off()

pdf(file = paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",
                  Sys.Date(),
                  "_",
                  phenotype,
                  "_",
                  enrich_fill,
                  ifelse(common_bg, "_common-bg", NA),
                  "_hm-enrich_motif-disruption.pdf"),
    width = 12,
    height = 11)
hm.enrich
dev.off()

ggsave(paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",
              Sys.Date(),
              "_",
              phenotype,
              "_",
              enrich_fill,
              ifelse(common_bg, "_common-bg", NA),
              "_dp-enrich_motif-disruption.pdf"),
       dp.enrich.egg,
       device = "pdf",
       width = 12,
       height = 11)


