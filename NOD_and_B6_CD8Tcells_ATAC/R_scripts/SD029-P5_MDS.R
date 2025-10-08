library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(dplyr)
library(egg)

# load data

proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
scratch.dir = "/scratch2/devoes/SD029"
half_width = 150
phenotypes = c("CART","LCMVcl13","Naive")

# metadata
meta = read.delim(paste0(proj.dir, "/metadata/low-noise_metadata.tsv"))

# limma voom model
lmv.v = readRDS(paste0(proj.dir, "/data/v_", half_width, "bp.rds"))


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

# Euclidean distance ####

# measures of similarity between samples: euclidean distance
lognorm.counts = lmv.v$E

colnames(lognorm.counts) = str_remove(colnames(lognorm.counts), pattern = "_REP.*")
sample_labels = paste(meta$Background, meta$Phenotype, meta$SampleID, sep = "_")
names(sample_labels) = meta$SampleID
colnames(lognorm.counts) = sample_labels[colnames(lognorm.counts)]

lognorm.counts = lognorm.counts[, grepl(paste0(phenotypes,collapse = "|"), colnames(lognorm.counts))]

meta = meta[grepl(paste0(phenotypes,collapse = "|"), meta$Phenotype),]
  
ecl.dist = t(lognorm.counts) %>% 
  dist(method = "euclidean") %>%
  as.matrix()

# Heatmap #### 

Heatmap(matrix = ecl.dist,
        col = hcl.colors(n = 100, palette = "YlOrBr", rev = F)[1:85], # YlOrBr, YlGr
        name = "Euclidean Distance",
        border = T,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8) ,
        heatmap_width = unit(6, "in"),
        heatmap_height = unit(6, "in"),
        #rect_gp = ,
        cluster_rows = T)
# row_dend_side = )


# MDS ####

# edgeR::plotMDS.DGEList() # would have been easier had I implemented this into my original lmv script

mds.fit = cmdscale(ecl.dist, k = 2, eig = T)
mds.df = mds.fit$points %>% as.data.frame() 
colnames(mds.df) = c("dim1", "dim2")
mds.df$Genotype =  mds.fit$points %>% rownames() %>% str_remove(pattern = "_.*")
mds.df$Phenotype = mds.fit$points %>% rownames() %>% 
  str_remove(pattern = "(B6|NOD)_") %>% str_remove(pattern = "_[0-9].*")
mds.df$Group = mds.fit$points %>% rownames() %>% str_remove(pattern = "_[0-9].*")

mds.colors = fill.colors[grep("B6", names(fill.colors))]
names(mds.colors) = names(mds.colors) %>% str_remove(pattern = "B6_")

p.mds =  ggplot(mds.df,
                aes(x = dim1, y = dim2, shape = Genotype, color = Phenotype)) + 
  geom_point(size = 2) +
  theme_classic() + 
  scale_color_manual(values = mds.colors) + 
  scale_shape_manual(values = c(15,0)) +
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

p.mds.egg = set_panel_size(p.mds, width = unit(2, "in"), height = unit(2,"in"))

# need to get legend to fit
# OR forgo margin and label with text in affinity?
ggsave(paste0(scratch.dir,"/figures/mds_cart_naive.pdf"),
       plot = p.mds.egg,
       width = 3,
       height = 3,
       device = "pdf")

ggsave(paste0(scratch.dir,"/figures/mds-wide_cart_naive.pdf"),
       plot = p.mds.egg,
       width = 3,
       height = 5,
       device = "pdf")


p.mds.group =  ggplot(mds.df,
                aes(x = dim1, y = dim2, color = Group)) + 
  geom_point(size = 2, shape = 0) +
  theme_classic() + 
  scale_color_manual(values = fill.colors) + 
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
p.mds.group

p.mds.group.egg = set_panel_size(p.mds.group, width = unit(2, "in"), height = unit(2,"in"))

# need to get legend to fit
# OR forgo margin and label with text in affinity?
ggsave(paste0(scratch.dir,"/figures/mds_cart_naive_group-square.pdf"),
       plot = p.mds.group.egg,
       width = 3,
       height = 3,
       device = "pdf")
