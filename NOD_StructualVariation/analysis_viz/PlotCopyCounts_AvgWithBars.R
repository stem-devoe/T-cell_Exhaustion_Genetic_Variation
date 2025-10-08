library(ggplot2)
library(dplyr)
library(stringr)
library(egg)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

# https://stackoverflow.com/questions/12399506/how-do-i-create-a-categorical-scatterplot-in-r-like-boxplots


metadata = read.delim("/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv")
plot_aes = read.delim("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/plot_aes.txt")
plot_shapes = setNames(plot_aes$Shape, plot_aes$Feature)
plot_colors = setNames(plot_aes$Color, plot_aes$peakID)

# denovo ####
denovo_libSize = read.csv("/Scottbrowne/members/smd/Projects/SD032/denovo_mapped_reads.csv", 
                          header = F)
denovo_libSize = denovo_libSize[!grepl("mapqgt1",denovo_libSize[,1]),]
denovo_libSize = data.frame(SampleID = basename(denovo_libSize[,1]) %>% str_remove(pattern = "_REP1.*") , 
                            libSize = denovo_libSize[,2])

copy_data = read.delim("/Scottbrowne/seq/tmp/devoes/SD032/CopyBiasCounts/copy_counts.tsv",
                       header = F)
# copy_data = read.delim("/Scottbrowne/seq/tmp/devoes/SD032/CopyBiasCounts_NonCNV_NoHetSNP/copy_counts.tsv",
#                        header = F)

colnames(copy_data) = c("copy_chr","copy_start","copy_end","peakID","bam","reads")

copy_data = copy_data %>% mutate(SampleID = basename(bam) %>% str_remove(pattern = "_REP1.*")) %>% left_join(copy_data, y= metadata, by = "SampleID") %>%
  left_join(y = denovo_libSize, by = "SampleID") %>% 
  mutate(cpm = reads/libSize*1e6) %>%
  mutate(copyID = paste(copy_chr,copy_start, sep = "_")) %>% 
  relocate(cpm, .after = reads) %>%
  left_join(y = plot_aes, by = "peakID")

copy_data.summary = copy_data %>%
  filter(Phenotype == "CART" & copy_chr == "OW971783.1") %>% 
  group_by(copy_start) %>% 
  summarize(mean_cpm = mean(cpm), 
            stdev = sd(cpm),
            se = plotrix::std.error(cpm), 
            ci.lower = t.test(cpm)$conf.int[1], 
            ci.upper = t.test(cpm)$conf.int[2])
copy_aes = copy_data[,c('copy_chr','copy_start','copyID','peakID','Feature','Color','Shape')] %>%
  filter(copy_chr == "OW971783.1") %>% 
  unique()
copy_data.summary = left_join(copy_data.summary, copy_aes, by = "copy_start")
plot_shapes['Cahan-CNV'] = 13
plot_colors_new = plot_colors
color_ind= grep(pattern = "#", plot_colors)
set.seed(101)
#plot_colors_new[color_ind] = rainbow(n = length(color_ind)) %>% sample()
#plot_colors_new[color_ind] = RColorBrewer::brewer.pal(n = length(color_ind), name = "Paired") %>% sample() # Dark2
plot_colors_new[color_ind] = hcl.colors(n = length(color_ind), palette = "Roma")  %>% sample() # Dynamic, Roma, 
plot_colors_new['chr1_171551434'] = "grey30"
p = copy_data.summary %>% 
  #filter(Feature != "Cahan-CNV") %>% 
  #filter(peakID == "chr1_171561902" | peakID == "chr1_171564324" ) %>%
  ggplot(aes(x = copy_start, y = mean_cpm, color = peakID, shape = Feature))+
  geom_point( size = 3) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper)) +
  #geom_errorbar(aes(ymin = mean_cpm - stdev, ymax = mean_cpm + stdev)) +
  #geom_errorbar(aes(ymin = mean_cpm - se, ymax = mean_cpm + se)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
  scale_shape_manual(values = plot_shapes) + 
  scale_color_manual(values = plot_colors_new) + 
 # xaxis.gap.fix +
  x_at_zero + 
  my.theme
p

p.egg = set_panel_size(p, width = unit(425, "mm"), height = unit(125,'mm')) #326.3, 425

ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_NOD_Cd244.pdf"),
       plot = p.egg,
       height = 10,
       width = 20,
       units = "in",
       device = "pdf")

# repeated measures anova (because mult measures within same sample) ####

# sampleID, copyID - variable, peakID - group, cpm - dependent Vari
copy_rma = copy_data[,c('copyID','peakID','SampleID','cpm','Phenotype')] %>% 
  filter(Phenotype == "CART" & (peakID == "chr1_171561902" | peakID == "chr1_171564324"))

# check assumptions
copy_rma %>%
  group_by(copyID, peakID) %>%
  summarize(mean = mean(cpm), sd = sd(cpm)) %>% 
  arrange(peakID) # do see heteroscedasticity so that means we should use the non parametric alternative

ggplot(copy_rma, aes(x = copyID, y = cpm, group = peakID)) +
  geom_point() + facet_grid(. ~ peakID) + theme_classic()

ggpubr::ggqqplot(copy_rma, "cpm", facet.by = c("peakID","copyID")) # deviation from normality on a couple but this is what violation test is most resistant too

copy_rma.l = split(copy_rma, f = as.factor(copy_rma$peakID))
  
# parametric (don't use)
# copy_rma_res = rstatix::anova_test(copy_rma.l[[1]], within = copyID, dv = cpm, wid = SampleID ) # sphericity is checked as part of this

# non-parametric
copy_friedman = lapply(copy_rma.l, FUN = rstatix::friedman_test, formula = cpm ~ copyID | SampleID)
copy_friedman

# mmarge ####
sd029_libSize = read.csv("/Scottbrowne/members/smd/Projects/SD032/mmarge_mapped_reads.csv", # MMARGE mapped and shifted
                         header = F)
sd029_libSize = data.frame(SampleID = basename(sd029_libSize[,1]) %>% str_remove(pattern = "_REP1.*") , 
                           libSize = sd029_libSize[,2])

sd029_counts = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/ArchR150bpExt_featureCounts.rds")
sd029_counts = sd029_counts$counts
sd029_counts = sd029_counts[unique(copy_data$peakID),] %>% t() %>% as.data.frame()
sd029_counts$bam = rownames(sd029_counts)
rownames(sd029_counts) = NULL
sd029_counts$SampleID = sd029_counts$bam %>% basename() %>% str_remove(pattern = "_REP1.*")
sd029_counts = reshape2::melt(sd029_counts, variable.name = "peakID", value.name = "reads")
sd029_counts = sd029_counts %>% # tidyr::pivot_longer(!SampleID, names_to = "peakID", values_to = "reads") %>% # pivot_longer was acting unexpectedly
  left_join( metadata, by = "SampleID") %>%
  left_join(y = sd029_libSize, by = "SampleID") %>% 
  mutate(cpm = reads/libSize*1e6)


sd029.summary = sd029_counts %>%
  filter(Phenotype == "CART") %>% 
  group_by(Strain, peakID) %>% 
  summarize(mean_cpm = mean(cpm), 
            stdev = sd(cpm),
            se = plotrix::std.error(cpm), 
            ci.lower = t.test(cpm)$conf.int[1], 
            ci.upper = t.test(cpm)$conf.int[2])


p.sd029 = sd029.summary %>% filter(peakID == "chr1_171561902" | peakID == "chr1_171564324") %>%
  ggplot(aes(Strain, mean_cpm, shape = Strain, color = Strain, group = peakID)) +
  geom_point( size = 3) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper)) +
  #geom_errorbar(aes(ymin = mean_cpm - stdev, ymax = mean_cpm + stdev)) +
  #geom_errorbar(aes(ymin = mean_cpm - se, ymax = mean_cpm + se)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ peakID, scales = "free") +
  scale_y_continuous(limits = c(0,max(copy_data.summary$ci.upper)), expand = expansion(mult = c(0,0.05))) + my.theme
p.sd029

p.sd029.egg = set_panel_size(p.sd029, width = unit(40,"mm"), height = unit(125, "mm"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_B6-NOD_Cd244_mm10.pdf"),
       plot = p.sd029.egg,
       height = 10,
       width = 8,
       device = "pdf")
