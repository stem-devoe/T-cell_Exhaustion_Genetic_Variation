# "/scratch2/devoes/SD029/hSNP"
library(stringr)
library(dplyr)
library(ggplot2)
library(egg)

# hsnp = read.delim("/scratch2/devoes/SD029/hSNP/mm10_1kb-window_200bp-step.hSNPcounts", header = F) # 5kb 1kb
# colnames(hsnp) = c("chr","start","stop","hSNP")
# quantile(hsnp$hSNP, probs = seq(0.96,1,0.001))
# 
# hsnp %>% filter(hSNP > 0) %>% # limit # of points to plot
#   mutate(rank = dense_rank(dplyr::desc(hSNP))) %>%
#   ggplot(aes(x = rank, y = hSNP)) + geom_point()


#hsnp2 = read.delim("/scratch2/devoes/SD029/hSNP/mm10_5kb-window_1kb-step.hSNPcounts", header = F) # 5kb 1kb
hsnp2 = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/mm10_5kb-window_1kb-step.hSNPcounts", header = F) # 5kb 1kb
colnames(hsnp2) = c("chr","start","stop","hSNP")
quantile(hsnp2$hSNP, probs = seq(0.9,1,0.005))
quantile(hsnp2$hSNP, probs = seq(0.9995,1,0.00001))
# hsnp2 %>%
#   mutate(rank = dense_rank(dplyr::desc(hSNP))) %>%
#   ggplot(aes(x = rank, y = hSNP)) + 
#   geom_point() +
#   geom_hline(yintercept = quantile(hsnp2$hSNP, probs = c(0.99998)), linetype = "dashed",color = "red") + 
#   geom_hline(yintercept = 320, color = "purple", linetype = "dashed")
# inflection at ~ 0.99998 (306 hSNP)
quantile(hsnp2$hSNP[hsnp2$hSNP > 0], probs = seq(0.9,1,0.01))

hdense = hsnp2[hsnp2$hSNP > quantile(hsnp2$hSNP[hsnp2$hSNP > 0], probs = 0.99) , 1:3]

options(scipen=9999)
write.table(hdense, 
            file = "/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/dense-hSNP_from_mm10_5kb-window_1kb-step.bed",
            #file = "/scratch2/devoes/SD029/hSNP/dense-hSNP_from_mm10_5kb-window_1kb-step.bed",
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

nonzero = hsnp2[hsnp2$hSNP > 0, ]
p.nz = nonzero %>%
  mutate(rank = dense_rank(dplyr::desc(hSNP))) %>%
  ggplot(aes(x = rank, y = hSNP)) +
  geom_point() +
  geom_hline(yintercept = quantile(nonzero$hSNP, probs = c(0.99)), linetype = "dashed",color = "red") + 
  theme_classic() + my.theme + xaxis.gap.fix
p.nz

p.nz.egg = set_panel_size(p.nz, width = unit(3, "in"), height = unit(3,"in"))
ggsave("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/nonzero_waterfall.pdf",
       plot = p.nz.egg,
       device = "pdf",
       unit = "in",
       width = 4,
       height = 4,
       dpi = 300)

head(hsnp2)
hsnp2.order = hsnp2[order(hsnp2$hSNP,decreasing = T),]
hsnp2.order$rank = seq_len(nrow(hsnp2.order))
head(hsnp2.order)

p.top = hsnp2.order[1:1e5,] %>%
  #hsnp2.order[1:round((nrow(hsnp2.order)/20)),] %>%
  ggplot(aes(x = rank, y = hSNP)) + 
  geom_point() +
  geom_hline(yintercept = quantile(nonzero$hSNP, probs = c(0.99)), linetype = "dashed",color = "red") + 
  theme_classic() + my.theme + xaxis.gap.fix
p.top
p.top.egg = set_panel_size(p.top, width = unit(3, "in"), height = unit(3,"in"))
ggsave("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/top100k_waterfall.pdf",
       plot = p.top.egg,
       device = "pdf",
       unit = "in",
       width = 4,
       height = 4,
       dpi = 300)

p.all = hsnp2.order %>%
  ggplot(aes(x = rank, y = hSNP)) +
  geom_point() +
  geom_hline(yintercept = quantile(nonzero$hSNP, probs = c(0.99)), linetype = "dashed",color = "red") + 
  theme_classic() + my.theme + xaxis.gap.fix
#p.all
p.all.egg = set_panel_size(p.all, width = unit(3, "in"), height = unit(3,"in"))
ggsave("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/all_waterfall.pdf",
       plot = p.all.egg,
       device = "pdf",
       unit = "in",
       width = 4,
       height = 4,
       dpi = 300)

# shell commands
# bedtools makewindows -g /Scottbrowne/members/smd/genomes/mm10/mm10.chrom.sizes -w 1000 -s 200 > /scratch2/devoes/SD029/hSNP/mm10_1kb-window_200bp-step.bed
# bedtools makewindows -g /Scottbrowne/members/smd/genomes/mm10/mm10.chrom.sizes -w 5000 -s 1000 > /scratch2/devoes/SD029/hSNP/mm10_5kb-window_1kb-step.bed
# bedtools intersect -c -a /scratch2/devoes/SD029/hSNP/mm10_1kb-window_200bp-step.bed -b /scratch2/devoes/SD029/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.hqHetSNP.vcf.bed > /scratch2/devoes/SD029/hSNP/mm10_1kb-window_200bp-step.hSNPcounts
# bedtools intersect -c -a /scratch2/devoes/SD029/hSNP/mm10_5kb-window_1kb-step.bed -b /scratch2/devoes/SD029/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.hqHetSNP.vcf.bed > /scratch2/devoes/SD029/hSNP/mm10_5kb-window_1kb-step.hSNPcounts