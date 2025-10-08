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

# denovo all peaks ####
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
  geom_point( size = 2, stroke = 0.5) +
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

width = 5
height = 2
p.egg = set_panel_size(p, width = unit(width, "in"), height = unit(height,'in'))

ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_",width,"x",height,"_THESIS_NOD_Cd244.pdf"),
       plot = p.egg,
       height = height + 4,
       width = width + 4,
       units = "in",
       device = "pdf")

# denovo select peaks all samples ####
p.select = copy_data %>%
  filter(copy_chr == "OW971783.1", Feature == "WithinCd244") %>% #, peakID == "chr1_171561902" | peakID ==  "chr1_171564324") %>% # | peakID ==  "chr1_171564324"
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
  facet_grid(rows = vars(Phenotype),cols = vars(peakID),scales = "free_x") #. ~ peakID*Phenotype
p.select
egg.len = 0.75
p.select.egg = set_panel_size(p.select, width = unit(egg.len, "in"), heigh = unit(egg.len, "in"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_",egg.len,"in_THESIS_Cd244-RegEle_byPeak.pdf"),
       plot = p.select.egg,
       height = 8,
       width = 10,
       units = "in",
       device = "pdf")

p.select.ylim = p.select + coord_cartesian(ylim = c(0,6.5))
p.select.ylim.egg = set_panel_size(p.select.ylim, width = unit(egg.len, "in"), heigh = unit(egg.len, "in"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_",egg.len,"in_THESIS_Cd244-RegEle_byPeak-SetYLims.pdf"),
       plot = p.select.ylim.egg,
       height = 8,
       width = 10,
       units = "in",
       device = "pdf")


# mmarge select peak all samples ####
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
  left_join(y = plot_aes, by = "peakID") %>%
  mutate(cpm = reads/libSize*1e6)


# sd029.summary = sd029_counts %>%
#   filter(Phenotype == "CART") %>% 
#   group_by(Strain, peakID) %>% 
#   summarize(mean_cpm = mean(cpm), 
#             stdev = sd(cpm),
#             se = plotrix::std.error(cpm), 
#             ci.lower = t.test(cpm)$conf.int[1], 
#             ci.upper = t.test(cpm)$conf.int[2])


p.sd029  = sd029_counts %>% 
  filter(Feature == "WithinCd244") %>%
  group_by(Phenotype, peakID) %>%
  ggplot(aes(x = Background, 
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
  facet_grid(rows = vars(Phenotype),cols = vars(peakID),scales = "free_x") +  #. ~ peakID*Phenotype
  coord_cartesian(ylim = c(0,6.5))

p.sd029.egg = set_panel_size(p.sd029, width = unit(egg.len,"in"), height = unit(egg.len, "in"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_",egg.len,"in_THESIS_Cd244-RegEle_byPeak_mm10.pdf"),
       plot = p.sd029.egg,
       height = 8,
       width = 10,
       device = "pdf")



# repeated measures anova (because mult measures within same sample) ####

# sampleID, copyID - variable, peakID - group, cpm - dependent Vari
copy_rma = copy_data[,c('copyID','peakID','SampleID','cpm','Phenotype','Feature')] %>% 
  filter(#Phenotype == "CART" & 
    Feature == "WithinCd244" # & (peakID == "chr1_171561902" | peakID == "chr1_171564324")
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

# parametric (don't use)
# copy_rma_res = rstatix::anova_test(copy_rma.l[[1]], within = copyID, dv = cpm, wid = SampleID ) # sphericity is checked as part of this

# non-parametric
copy_friedman = lapply(copy_rma.l, FUN = rstatix::friedman_test, formula = cpm ~ copyID | SampleID)
copy_friedman
lapply(copy_friedman, FUN = function(x){x$p}) %>% p.adjust(method = "fdr")



# bar graph num blast hits #####

# need to summarize counts per nod and per b6 for peaks we ran blast on
num_copies = read.csv("/Scottbrowne/seq/tmp/devoes/SD032/blast_all_peaks/hits_per_peak.csv", header = F)
colnames(num_copies) = c("peakID","Genome","NumBlastHits")
num_copies = num_copies %>% 
  filter(peakID %in% c(plot_aes$peakID,
                       "chr1_171606632","chr1_171606972","chr1_171607316","chr1_171610071","chr1_171610814")) %>% 
  left_join(y = plot_aes, by = "peakID")

p.nblast = num_copies %>%
  mutate(Feature = factor(Feature, c("Upstream","Cahan-CNV","WithinCd244","Downstream"))) %>%
 # filter() %>%
  ggplot(aes(fill = Genome, x = peakID, y = NumBlastHits)) +
  geom_bar(position = "dodge",stat = "identity") +
 # facet_grid(. ~ Feature, scales = "free_x") +
  theme_classic() +
  scale_fill_manual(values = c("black","grey50")) +
  theme(axis.text.x = element_text(angle = 90)) +
  my.theme +
  xaxis.gap.fix
p.nblast.egg = set_panel_size(p.nblast, width = unit(3,"in"), height = unit(2,"in"))
ggsave(
  paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_THESIS_Cd244-RegEle_NumBlastHits_noFacet.pdf"),
  plot = p.nblast.egg,
  height = 6,
  width = 6,
  device = "pdf"
)


# imputed vs de novo ####
# get max ...902 and 324 NOD copy for CART because of num samples


# identify copy with highest mean_cpm 
# (calc of mean_cpm was already done and filtered for CART and hits on std chr)
max.copy.summary = copy_data.summary %>%
  filter(peakID == "chr1_171561902" | peakID == "chr1_171564324") %>% 
  group_by(peakID) %>%
  filter(mean_cpm == max(mean_cpm))

max.copy = copy_data %>% filter(copyID %in% max.copy.summary$copyID & Phenotype == "CART")

# get mean values
max.copy.summary
sd029_counts %>% filter(Phenotype == "CART" & (peakID == "chr1_171561902" | peakID == "chr1_171564324")) %>%
  group_by(peakID,Background) %>% 
  summarize(mean_cpm = mean(cpm),
            stdev = sd(cpm),
            se = plotrix::std.error(cpm), 
            ci.lower = t.test(cpm)$conf.int[1],
            ci.upper = t.test(cpm)$conf.int[2]) %>% print(n=Inf)


# check dist
ggplot(max.copy, aes(x = copy_start, y = cpm, group = peakID)) + geom_point() + geom_boxplot()
max.copy %>%  ggpubr::ggqqplot("cpm", facet.by = "peakID") # 902 deviates a bit from normal
max.copy %>% mutate(cpm = log2(cpm + 1e-7)) %>% ggpubr::ggqqplot("cpm", facet.by = "peakID")

sd029_counts %>%  filter(peakID == "chr1_171561902" | peakID == "chr1_171564324") %>% 
  filter(Phenotype == "CART") %>%
  ggpubr::ggqqplot("cpm", facet.by = c("peakID","Background"))
  #ggplot(aes(x = peakID, y = cpm, group = peakID)) + geom_point() + geom_boxplot()


## man whitney U ####

# peakID, genome, cpm
artefact.df = sd029_counts[,c('SampleID','cpm','peakID','Phenotype','Background')] %>%  filter(peakID == "chr1_171561902" | peakID == "chr1_171564324") %>% 
  filter(Phenotype == "CART") 
artefact.df$copyID = artefact.df$peakID
artefact.df$Genome = artefact.df$Background %>% str_replace(pattern = "NOD",replacement = "NOD_MMARGE")
max.copy$Genome = "NOD_SHILTJ_v3"
artefact.df = rbind(artefact.df, max.copy[,c('SampleID','cpm','peakID','Phenotype','Background','copyID','Genome')])

runMWU = function(df){
  mmarge = df[df$Genome == "NOD_MMARGE",]
  b6 =  df[df$Genome == "B6",]
  denovo = df[df$Genome == "NOD_SHILTJ_v3",]
  
  nodvmmarge= wilcox.test(denovo$cpm, mmarge$cpm,  alternative = "two.sided")
  nodvmm10 = wilcox.test(denovo$cpm, b6$cpm, alternative = "two.sided")
  mmargevmm10 = wilcox.test(mmarge$cpm, b6$cpm, alternative = "two.sided")
  res = list(nodvmmarge, nodvmm10, mmargevmm10)
  names(res) = c("DenovoVsMMARGE","DenovoVsB6","MMargeVsB6")
  return(res)
}
artefact.df %>% split( f = ~ peakID) %>% lapply( FUN = runMWU)

## plot ####

poi = c("chr1_171561902","chr1_171564324")
coi = max.copy$copyID %>% unique()
p.art = artefact.df %>% 
  group_by(peakID) %>%
  ggplot(aes(x = paste0(Genome,":",copyID) %>%
               factor(levels = c(paste0("B6",":",poi),
                                 paste0("NOD_SHILTJ_v3",":",coi),
                                 paste0("NOD_MMARGE",":",poi))), 
             y = cpm, 
             color = peakID, 
            # group = Phenotype, 
             shape = Background))+
  geom_point(position=position_jitter(width=0.1, height=0), 
             size = 1, 
             #shape = 0, sape by genome
             stroke = 0.5) +
  stat_summary(geom = "errorbar", width = 0.5,
               color = "black", linewidth = 0.5) +
  stat_summary(fun=mean,fun.max = mean, fun.min = mean, 
               geom="errorbar", width = 0.25, 
               color = "black", linewidth = 0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = plot_colors_new) + 
  # xaxis.gap.fix +
  x_at_zero + 
  my.theme + 
  facet_grid(.~peakID,scales = "free_x") +  #. ~ peakID*Phenotype
  coord_cartesian(ylim = c(0,6.5))

p.art.egg = set_panel_size(p.art, width = unit(egg.len,"in"), height = unit(egg.len,"in"))
ggsave(paste0("/Scottbrowne/members/smd/Projects/SD032/CopyBiasCounts/",Sys.Date(),"_",egg.len,"in_THESIS_Cd244-RegEle_artefact.pdf"),
       plot = p.art.egg,
       height = 4,
       width = 4,
       device = "pdf")
