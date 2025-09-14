library(ggplot2)
library(stringr)
library(dplyr)
library(egg)
library(GRaNIE)

peaks = c(
  "chr3_37387033",
  "chr3_37273703",
  "chr1_171502083",
  "chr1_171561868",
  "chr1_128611823",
  "chr4_136089879",
  "chr9_118519148",
  "chr4_150915264"
          )
#peaks = paste0("^",peaks,"$")

proj_path = "/Scottbrowne/members/smd/Projects/SD030"
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                     "2024-09-01","_peak-gene_pearson_GRN.rds"))
group = with(GRN@data$metadata, paste(Publication,TCR,Phenotype,Setting,Treatment4GroupLabel, sep = "_")) # peaks were called on proper groups, the meta data for peak group is wrong only in granie object and metadata file

labels = read.delim("/Scottbrowne/members/smd/Projects/SD030/sample_inputs/Group2Label.txt",
                    header = F)
labels = setNames(labels$V2, labels$V1)
sort_key = as.vector(labels)

group_label = labels[group]

atac = getCounts(GRN, type = "peaks")
rownames(atac) = atac$peakID %>% str_remove(pattern = "-.*") %>% 
  str_replace(pattern = ":", replacement = "_")

atac = atac[peaks ,
          !(grepl("peakID",colnames(atac)))] 
atac = atac %>% t() %>% as.data.frame()
atac$group_label = group_label

atac_mean = atac %>% group_by(group_label) %>% summarise(across(.cols = everything(), mean))


atac_mean_mtx = as.matrix(atac_mean[,!grepl("group",colnames(atac_mean))])
rownames(atac_mean_mtx) = atac_mean$group_label
sort_order = lapply(sort_key, 
                    FUN = function(key,groups){grep(key,groups, 
                                                    ignore.case = T, 
                                                    value = T)},
                      rownames(atac_mean_mtx)) %>% unlist()
atac_mean_mtx = atac_mean_mtx[sort_order,]
atac_mean_mtx_center = scale(atac_mean_mtx, center = T, scale = F)


# dot bar graph ####

atac.dot.df = atac_mean_mtx_center %>% as.data.frame() %>% mutate(group = rownames(atac_mean_mtx_center))
atac.dot.df = reshape2::melt(atac.dot.df)
atac.dot.df$group = factor(atac.dot.df$group, levels = rev(sort_order))

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

p.bardot = ggplot(atac.dot.df, aes(x = value, y = group, group = variable)) + 
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey40") + 
  geom_point(size = 2, col = "black") + facet_grid(. ~ variable, scales = "free") +
  theme_classic() + my.theme + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) # expand ranges so dots aren't sitting on edge of polot 
p.bardot

p.egg = set_panel_size(p.bardot, width = unit(45, "mm"), height = unit(183.7, "mm"))

ggsave(paste0(proj_path,"/figures/atac-signal.pdf"),
       plot = p.egg,
       width = 24,
       height = 12,
       device = "pdf")
