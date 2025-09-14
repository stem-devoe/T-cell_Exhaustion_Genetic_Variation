library(ggplot2)
library(stringr)
library(dplyr)
library(egg)
library(GRaNIE)

proj_path = "/Scottbrowne/members/smd/Projects/SD030"
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                    "2024-10-06","_peak-gene_pearson_GRN.rds"))
group = GRN@data$metadata$archr_peak_group # peaks were called on proper groups, the meta data for peak group is wrong only in granie object and metadata file

labels = read.delim("/Scottbrowne/members/smd/Projects/SD030/sample_inputs/Group2Label.txt",
                    header = F)
labels = setNames(labels$V2, labels$V1)
  
mds_label = labels[group]

study = GRN@data$metadata$Publication

rna = getCounts(GRN, type = "rna")
rownames(rna) = rna$ENSEMBL
rna = rna[,!grepl("ENSEMBL", colnames(rna))]
atac = getCounts(GRN, type =  "peaks")
rownames(atac) = atac$peakID
atac = atac[,!grepl("peakID",colnames(atac))]

# MDS RNA ####

ecl.dist = t(rna) %>% # samples as rows 
  dist(method = "euclidean") %>%
  as.matrix()

mds.fit = cmdscale(ecl.dist, k = 2, eig = T)
mds.df = mds.fit$points %>% as.data.frame() 
colnames(mds.df) = c("dim1", "dim2")

# fix metadata for SD030 samples: shape for publication, color for rest (or vice versa?)
mds.df = data.frame(mds.df, group, study, mds_label)
  
p.mds =  ggplot(mds.df,
                aes(x = dim1, 
                    y = dim2, 
                    shape = study, 
                    color = mds_label)) + 
  geom_point(size = 2) + 
  theme_classic() + 
 # scale_color_manual(values = mds.colors) + 
 # scale_shape_manual(values = c(15,0)) +
  theme(panel.border = element_rect(linewidth = 0.2, fill = NA),
        axis.line = element_line(linewidth = 0.2),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 7),
        #axis.ticks = element_line(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,1,0),
        legend.spacing = unit(0, "mm"), # spacing between legends
        #legend.key.height = ,
        legend.key.width = unit(0, "in"),
        legend.box.spacing = unit(0, "mm"),
        legend.title = element_blank()
  )

p.mds

p.mds.egg = set_panel_size(p.mds, width = unit(2, "in"), height = unit(2,"in"))

# need to get legend to fit
# OR forgo margin and label with text in affinity?
ggsave(paste0(proj_path,"/figures/rna_mds.pdf"),
       plot = p.mds.egg,
       width = 5,
       height = 5,
       device = "pdf")

ggsave(paste0(proj_path,"/figures/rna_mds-wide.pdf"),
       plot = p.mds.egg,
       width = 14,
       height = 5,
       device = "pdf")


# MDS ATAC ####

ecl.dist.atac = t(atac) %>% # samples as rows 
  dist(method = "euclidean") %>%
  as.matrix()

mds.atac.fit = cmdscale(ecl.dist.atac, k = 2, eig = T)
mds.atac.df = mds.atac.fit$points %>% as.data.frame() 
colnames(mds.atac.df) = c("dim1", "dim2")

mds.atac.df = data.frame(mds.atac.df, group, study, mds_label)

p.mds.atac =  ggplot(mds.atac.df,
                     aes(x = dim1, 
                         y = dim2, 
                         shape = study, 
                         color = mds_label)) + 
  geom_point(size = 2) + 
  theme_classic() + 
  theme(panel.border = element_rect(linewidth = 0.2, fill = NA),
        axis.line = element_line(linewidth = 0.2),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        title = element_text(size = 7),
        #axis.ticks = element_line(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,1,0),
        legend.spacing = unit(0, "mm"), # spacing between legends
        #legend.key.height = ,
        legend.key.width = unit(0, "in"),
        legend.box.spacing = unit(0, "mm"),
        legend.title = element_blank()
  )

p.mds.atac

p.mds.atac.egg = set_panel_size(p.mds.atac, width = unit(2, "in"), height = unit(2,"in"))

# need to get legend to fit
# OR forgo margin and label with text in affinity?
ggsave(paste0(proj_path,"/figures/atac_mds.pdf"),
       plot = p.mds.atac.egg,
       width = 5,
       height = 5,
       device = "pdf")

ggsave(paste0(proj_path,"/figures/atac_mds-wide.pdf"),
       plot = p.mds.atac.egg,
       width = 14,
       height = 5,
       device = "pdf")





