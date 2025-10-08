library(dplyr)
library(stringr)
library(ggplot2)
#library(UpSetR) # ComplexHeatmap is easier to use
library(ComplexHeatmap) # https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
library(egg) # fixes panel size for ggplots

# load data

fold_change = 1.5
half_width = 150
padj = 0.1
proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
scratch.dir = "/scratch2/devoes/SD029"

# metadata
genome = "mm10"
phenos = c("CART","LCMVcl13","Naive")

# peaks
peak.consensus.file = paste0(scratch.dir,
                             "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt",
                             "/globallambda/150bp",
                             "/All_Samples.fwp.filter.non_overlapping.bed")
peak.consensus = readRDS(paste0(scratch.dir,
                                "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt",
                                "/globallambda/150bp",
                                "/All_Samples.fwp.filter.non_overlapping.rds"))
peak.sig.files = list.files(paste0(proj.dir, "/data/",half_width,"bp/padj_",padj,"/sig_peaks/FC_", fold_change),
                            pattern = "up.bed|down.bed",
                            full.names = T) %>% grep(pattern = paste0(phenos, collapse = "|"), value = T) %>%
  grep(pattern = "NOD-B6", value = T)
peak.sig = lapply(peak.sig.files, FUN = rtracklayer::import)
names(peak.sig) = basename(peak.sig.files) %>% str_remove(pattern = ".bed")


# color palettes 
sort.key = c("B6_Naive","NOD_Naive",
             "B6_LCMVcl13","NOD_LCMVcl13",
             "B6_CART","NOD_CART")

fill.colors = c(B6_Naive = "#ff33ff",
                NOD_Naive = "#cc99ff",
                B6_LCMVcl13 = "#009900",
                NOD_LCMVcl13 = "#99cc66",
                #B6_CART = "#3333ff",
                B6_CART = "#4969ff",
                NOD_CART = "#6699ff")

# Upset plots ####

# get list of peak IDs per comparison
peakIDs.comparison <- lapply(peak.sig, FUN = function(x){x$name})

# rename list
names(peakIDs.comparison) <- names(peakIDs.comparison) %>% 
  str_remove(pattern = "NOD-B6_") %>% 
  str_replace(pattern = "down", replacement = " < ") %>%
  str_replace(pattern = "up", replacement = " > ") %>% 
  str_remove(pattern= "\\.") %>%
  paste0("NOD ", ., "B6 ", str_remove(string = ., pattern = "NOD | < | > | B6"))

# convert to binary matrix
peaks.mtx = list_to_matrix(peakIDs.comparison)

# make combination matrix
comb.mtx = make_comb_mat(peakIDs.comparison, mode = "distinct")

# plot
p.upset <- UpSet(comb.mtx,
                 heatmap_height = unit(2.5, "in"),
                 heatmap_width = unit(5, "in"),
                 pt_size = unit(3, "mm"),
                 lwd = 1,
                 set_order = seq_along(set_size(comb.mtx)), # I want the order I manually had
                 comb_order = order(comb_size(comb.mtx), decreasing = T), # do we order that makes sense by connected dots or by most to least 
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize =7,
                                     col = fill.colors[c('B6_CART','NOD_CART',
                                                         "B6_LCMVcl13","NOD_LCMVcl13",
                                                         'B6_Naive','NOD_Naive')]%>%as.vector()),
                 top_annotation = upset_top_annotation(comb.mtx, 
                                                       annotation_name_gp = gpar(fontsize = 0),
                                                       annotation_name_rot = 90,
                                                       add_numbers = T,
                                                       axis = F,
                                                       numbers_gp = gpar(fontsize = 6),
                                                       numbers_rot = 90),
                 right_annotation = upset_right_annotation(comb.mtx,
                                                           annotation_name_gp = gpar(fontsize = 0),
                                                           gp = gpar(fill = fill.colors[c('B6_CART','NOD_CART',
                                                                                          "B6_LCMVcl13","NOD_LCMVcl13",
                                                                                          'B6_Naive','NOD_Naive')],
                                                                     col = fill.colors[c('B6_CART','NOD_CART',
                                                                                         "B6_LCMVcl13","NOD_LCMVcl13",
                                                                                         'B6_Naive','NOD_Naive')]),
                                                           add_numbers = T,
                                                           axis = F,
                                                           numbers_gp = gpar(fontsize = 6),
                                                           bar_width = 0.5) 
)

p.upset

pdf(file = "/scratch2/devoes/SD029/figures/distinct-upset.pdf",
    width = 5.5,
    height = 3)
p.upset
dev.off()
