library(ggplot2)
library(dplyr)
library(stringr)
library(egg)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

# https://stackoverflow.com/questions/12399506/how-do-i-create-a-categorical-scatterplot-in-r-like-boxplots

avg_phenotype = "CART"
metadata = read.delim("/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv")
plot_aes = read.delim("/Scottbrowne/members/smd/Projects/SD032/Ccl27_Ccl21_Il11ra/plot_aes.txt")
plot_shapes = setNames(plot_aes$Shape, plot_aes$Feature)
plot_colors = setNames(plot_aes$Color, plot_aes$peakID)
a_peaks = plot_aes$peakID

# denovo all peaks ####
denovo_libSize = read.csv("/Scottbrowne/members/smd/Projects/SD032/denovo_mapped_reads.csv", 
                          header = F)
denovo_libSize = denovo_libSize[!grepl("mapqgt1",denovo_libSize[,1]),]
denovo_libSize = data.frame(SampleID = basename(denovo_libSize[,1]) %>% str_remove(pattern = "_REP1.*") , 
                            libSize = denovo_libSize[,2])

copy_data = read.delim("/Scottbrowne/seq/tmp/devoes/SD032/Ccl27_Ccl21_Il11ra/peaks_blast/NOD/copy_counts.tsv",
                       header = F)

colnames(copy_data) = c("copy_chr","copy_start","copy_end","peakID","strand","hitNum","bam","reads")

copy_data = copy_data %>% mutate(SampleID = basename(bam) %>% str_remove(pattern = "_REP1.*")) %>% left_join(copy_data, y= metadata, by = "SampleID") %>%
  left_join(y = denovo_libSize, by = "SampleID") %>% 
  mutate(cpm = reads/libSize*1e6) %>%
  mutate(copyID = paste(copy_chr, copy_start, sep = "_")) %>% 
  relocate(cpm, .after = reads) %>%
  left_join(y = plot_aes, by = "peakID")

copy_data.summary = copy_data %>%
  filter(Phenotype == avg_phenotype & copy_chr == "OW971786.1") %>% 
  group_by(copy_start) %>% 
  summarize(mean_cpm = mean(cpm), 
            stdev = sd(cpm),
            se = plotrix::std.error(cpm), 
            ci.lower = t.test(cpm)$conf.int[1], 
            ci.upper = t.test(cpm)$conf.int[2])
copy_aes = copy_data[,c('copy_chr','copy_start','copyID','peakID','Feature','Color','Shape')] %>%
  filter(copy_chr == "OW971786.1") %>% 
  unique()
copy_data.summary = left_join(copy_data.summary, copy_aes, by = "copy_start")
plot_colors_new = plot_colors
color_ind= grep(pattern = "#", plot_colors)
set.seed(101)
#plot_colors_new[color_ind] = rainbow(n = length(color_ind)) %>% sample()
#plot_colors_new[color_ind] = RColorBrewer::brewer.pal(n = length(color_ind), name = "Paired") %>% sample() # Dark2
plot_colors_new[color_ind] = hcl.colors(n = length(color_ind), palette = "Roma")  %>% sample() # Dynamic, Roma, 

p = copy_data.summary %>% 
  filter(peakID %in% a_peaks) %>% 
  ggplot(aes(x = copy_start, y = mean_cpm, color = peakID, shape = Feature))+
  geom_point(size = 2, stroke = 0.5) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), linewidth = 0.5) +
  #geom_errorbar(aes(ymin = mean_cpm - stdev, ymax = mean_cpm + stdev)) +
  #geom_errorbar(aes(ymin = mean_cpm - se, ymax = mean_cpm + se)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
  scale_shape_manual(values = plot_shapes) + 
  scale_color_manual(values = plot_colors_new) + 
  # xaxis.gap.fix +
  x_at_zero + 
  my.theme 
p

width = 4
height = 2
p.egg = set_panel_size(p, width = unit(width, "in"), height = unit(height,'in'))

ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/Ccl27_Ccl21_Il11ra/",
              Sys.Date(),"_",width,"x",height,"_",avg_phenotype,"_NOD_THESIS_Il11ra.pdf"),
       plot = p.egg,
       height = height + 4,
       width = width + 4,
       units = "in",
       device = "pdf")

# denovo select peaks all samples ####
p.select = copy_data %>%
  filter(copy_chr == "OW971786.1", peakID %in% a_peaks) %>% 
  group_by(copy_start, Phenotype, peakID) %>%
  ggplot(aes(x = copyID, 
             y = cpm, 
             color = peakID, 
             group = Phenotype, 
             shape = Feature))+
  geom_point(position=position_jitter(width=0.1, height=0), size = 1, shape = 0, stroke = 0.5) +
  # geom_errorbar(stat='summary', width= 0.5, color = "black") + # mean and SE
  stat_summary(geom = "errorbar", width = 0.5,
               color = "black", linewidth = 0.5) +
  stat_summary(fun=mean,fun.max = mean, fun.min = mean, 
               geom="errorbar", width = 0.25, 
               color = "black", linewidth = 0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90)) +
  # scale_shape_manual(values = plot_shapes) + 
  scale_color_manual(values = plot_colors_new) + 
  # xaxis.gap.fix +
  x_at_zero + 
  my.theme + 
  facet_grid(
    rows = vars(Phenotype),cols = vars(peakID),
    # rows = vars(peakID), cols = vars(Phenotype),
    space = "free_x",
    scales = "free") #. ~ peakID*Phenotype
p.select
egg.len = 0.75
p.select.egg = set_panel_size(p.select, width = unit(egg.len, "in"), heigh = unit(egg.len, "in"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/Ccl27_Ccl21_Il11ra/",
              Sys.Date(),"_",egg.len,"in_NOD_THESIS_Il11ra-RegEle_byPeak.pdf"),
       plot = p.select.egg,
       height = 8,
       width = 10,
       units = "in",
       device = "pdf")

# p.select.ylim = p.select + coord_cartesian(ylim = c(0,6.5))
# p.select.ylim.egg = set_panel_size(p.select.ylim, width = unit(egg.len, "in"), heigh = unit(egg.len, "in"))
# ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/Ccl27_Ccl21_Il11ra/",Sys.Date(),"_",egg.len,"in_NOD_THESIS_Il11ra-RegEle_byPeak-SetYLims.pdf"),
#        plot = p.select.ylim.egg,
#        height = 8,
#        width = 10,
#        units = "in",
#        device = "pdf")



# repeated measures anova (because mult measures within same sample) ####

# sampleID, copyID - variable, peakID - group, cpm - dependent Vari
copy_rma = copy_data[,c('copyID','peakID','SampleID','cpm','Phenotype','Feature')] %>% 
  filter(#Phenotype == "CART" & 
    peakID %in% a_peaks &
      Feature != "Outside"
  )

# check assumptions
copy_rma %>%
  group_by(copyID, peakID,Phenotype) %>%
  summarize(mean = mean(cpm), sd = sd(cpm)) %>% 
  arrange(peakID) # do see heteroscedasticity so that means we should use the non parametric alternative

ggplot(copy_rma, aes(x = copyID, y = cpm, group = peakID)) +
  geom_point() + facet_grid(rows = vars(Phenotype), cols = vars(peakID)) + theme_classic()

copy_rma %>% filter(Phenotype == "CART") %>%
  ggpubr::ggqqplot( "cpm", facet.by = c("peakID","copyID")) # deviation from normality on a couple but this is what violation test is most resistant too
copy_rma %>% filter(Phenotype == "LCMVcl13") %>%
  ggpubr::ggqqplot( "cpm", facet.by = c("peakID","copyID")) # deviation from normality on a couple but this is what violation test is most resistant too

copy_rma.l = split(copy_rma, f = as.factor(paste0(copy_rma$peakID,"_", copy_rma$Phenotype)))


# non-parametric
copy_friedman = lapply(copy_rma.l, FUN = rstatix::friedman_test, formula = cpm ~ copyID | SampleID)
copy_friedman
lapply(copy_friedman, FUN = function(x){x$p}) %>% p.adjust(method = "fdr")

