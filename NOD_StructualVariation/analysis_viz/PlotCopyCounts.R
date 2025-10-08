library(ggplot2)
library(dplyr)
library(stringr)
library(egg)

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

copy_data %>% filter(Phenotype == "CART" & copy_chr == "OW971783.1" & Feature != "Cahan-CNV") %>% 
  ggplot( aes(x = copy_start, y = cpm, color = peakID, shape = Feature))+
  geom_jitter(position=position_jitter(width=0.1), size = 2) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
  scale_shape_manual(values = plot_shapes) + scale_color_manual(values = plot_colors)


## 
p.902 = copy_data %>% filter(Phenotype == "CART" & peakID == "chr1_171561902") %>% # split_by?
  ggplot( aes(x = copyID, y = cpm, shape = SampleID)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90)) +
  scale_shape_manual(values = c(15,16,17,0,1,2,3,8)) # 13

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_chr1_171561902_copies.pdf"),
       plot = p.902,
       height = 5,
       width = 7,
       units = "in",
       device = "pdf")

p.324 = copy_data %>% filter(Phenotype == "CART" & peakID == "chr1_171564324") %>% # split_by?
  ggplot( aes(x = copyID, y = cpm, shape = SampleID)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90)) +
  scale_shape_manual(values = c(15,16,17,0,1,2,3,8)) # 13

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_chr1_171564324_copies.pdf"),
       plot = p.324,
       height = 5,
       width = 7,
       units = "in",
       device = "pdf")

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

## It might be better to take normalized counts from v$E to display the NODvsB6 consensus peak because the cpm doesn't really capture it

sd029_counts %>% filter(Phenotype == "CART" & peakID == "chr1_171561902") %>%
  ggplot( aes(Strain, cpm, shape = Strain, color = Strain)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90))

sd029_counts %>% filter(Phenotype == "CART" & peakID == "chr1_171564324") %>%
  ggplot( aes(Strain, cpm, shape = Strain, color = Strain)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90))

voom = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/v_150bp.rds")
sd029_lmv = voom$E
sd029_lmv = sd029_lmv[unique(copy_data$peakID),] %>% t() %>% as.data.frame()
sd029_lmv$bam = rownames(sd029_lmv)
rownames(sd029_lmv) = NULL
sd029_lmv$SampleID = sd029_lmv$bam %>% basename() %>% str_remove(pattern = "_REP1.*")
sd029_lmv = reshape2::melt(sd029_lmv, variable.name = "peakID", value.name = "lmv_norm") %>% 
  left_join( metadata, by = "SampleID")

sd029_lmv %>% filter(Phenotype == "CART" & peakID == "chr1_171561902") %>%
  ggplot( aes(Strain, lmv_norm, shape = Strain, color = Strain)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90)) +
  scale_y_continuous(expand = c(0, 1.1), limits = c(0, NA)) +
  scale_color_manual(values = c("black","grey50"))

sd029_lmv %>% filter(Phenotype == "CART" & peakID == "chr1_171564324") %>%
  ggplot( aes(Strain, lmv_norm, shape = Strain, color = Strain)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  theme_classic() + theme(axis.text.x = element_text(angle =90)) +
  scale_y_continuous(expand = c(0, 1.1), limits = c(0, NA))


