library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
# updated 20240903 to use counts filtered decide test from phenotype comparisons (within same strain)

motif_db = "HOCOMOCO"
phenotype = "Naive" # CART LCMV Naive
save_date = "2024-09-04"
motif_family = read.delim("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_tf-info.txt",
                          header = T)
colnames(motif_family)[1] = "motif"
lfc.res = readRDS(file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                                motif_db, 
                                "/150bp/padj_0.1/FC_1.5/",
                                phenotype,
                                "/", 
                                save_date,
                                "_NoHetSNP_lfc-res_2sample.rds")
)

label_bool =  (lfc.res$p.adj < 0.05) & (abs(lfc.res$diff_mean_LFC) > 0.5)
lfc.res[label_bool,'motif']
lfc.res[lfc.res$p.adj < 0.05,]

lfc.res = left_join(lfc.res, motif_family, by = "motif")

lfc.res$Plot_Motif_Color[grep(pattern = "ETS|Runt|bHLH|bZIP|HMG|NR|RHR|T-box|Zf",lfc.res$Plot_Motif_Color,invert = T, ignore.case = T)] = NA

lfc.res$Point_Label = rep(NA, nrow(lfc.res))
lfc.res$Point_Label[label_bool] = lfc.res$motif[label_bool] %>% str_remove(pattern = "_MOUSE.*")

lfc.res[which(lfc.res$p.adj < 0.05 & lfc.res$diff_mean_LFC < -0.5),c(1:6,8)]

# scatter plot
# x = shared
# y = B6
# size = pval
# color = motif family

ggplot(lfc.res, aes(x = mean_shared_lfc, y = mean_B6_lfc, size = -log10(p.adj), color = Plot_Motif_Color)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label), size = 3, min.segment.length = 0, color = "black") +
  theme_classic()

# volcano plot
# x = diff LFC
# y = -log10padj
# color = motif family
pt.size = 0.3 
axis.size = 9 #7
label.size = 9 #6
title.size = 9 # 7
border.width = 0.2

p.volc = ggplot(lfc.res, aes(x = diff_mean_LFC, y = -log10(p.adj), color = Plot_Motif_Color, label = Point_Label)) + 
  geom_point(size = 1 , shape = 16) +
  geom_text_repel(size = 2, 
                  min.segment.length = 0 ,
                  segment.size = unit(0.05, "cm"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = Inf)
  ) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red', linewidth = border.width) + 
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'red', linewidth = border.width) + 
  geom_vline(xintercept = -0.5, linetype = 'dashed', color = 'red', linewidth = border.width) +
  xlab("difference in mean LFC") +
  theme(panel.border = element_rect(colour = "black", fill = NA, # add border to panel
                                    linewidth = border.width),
        legend.position = "none",
        axis.line = element_line(linewidth = border.width),
        axis.title = element_text(size = axis.size, margin = margin(0,0,0,0)), # axis title
        axis.text = element_text(size = label.size, margin = margin(0,0,0,0)), # tickmark labels
        plot.title = element_text(size = title.size, margin = margin(0,0,0,0)), # plot title
        axis.ticks.length = unit(0.05, "cm"), # adjust length of tick (outside panel)
        axis.ticks =  element_line(linewidth = border.width)) # adjust width of tick(outside panel) (mm)

p.volc
# orig 3x3
p.2 = egg::set_panel_size(p.volc, width = unit(2, "in"), height = unit(2,"in"))
p.3 = egg::set_panel_size(p.volc, width = unit(3, "in"), height = unit(3,"in"))

ggsave(paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",Sys.Date(),"_", phenotype, "_3in_NoHetSNP_motifdisruption_volcano.pdf"),
       plot = p.3,
       device = "pdf",
       unit = "in",
       width = 4,
       height = 4)

ggsave(paste0("/Scottbrowne/seq/tmp/devoes/SD029/figures/",Sys.Date(),"_", phenotype, "_2in_NoHetSNP_motifdisruption_volcano.pdf"),
       plot = p.2,
       device = "pdf",
       unit = "in",
       width = 4,
       height = 4)


# get values for venn diagram ####
half_width=150
padj=0.1
fold_change=1.5
proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
rm_peak = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/ConsensusPeaks_intersect_AnyHetSNP.bed",
                     header = F)
rm_peak = rm_peak$V4
#dt.res = readRDS(paste0(proj.dir, "/data/",half_width, "bp/padj_",padj,"/decideTests-Phenotype_FC",fold_change,".rds"))
dt.res = readRDS(paste0(proj.dir, "/data/",half_width, "bp/padj_",padj,"/decideTests-Phenotype_FC",fold_change,"_countsFiltered.rds"))
dt.res = dt.res[!grepl(pattern = paste0(rm_peak,collapse = "|"),rownames(dt.res)),]

if(phenotype != "Naive"){
  shared.pheno = limma::vennCounts(dt.res[,grepl(phenotype, colnames(dt.res))], include = "up")
  shared.pheno.df = data.frame(change = c("B6","NOD","Shared"), 
                               num_peaks = shared.pheno[2:4,3])
  print(phenotype)
  shared.pheno.df
} else{
  shared.pheno = limma::vennCounts(dt.res[,grepl(phenotype, colnames(dt.res))], include = "down")
  print(phenotype)
  shared.pheno
}
