#stat_ecdf
#ecdf

library(ggplot2)
library(dplyr)
library(stringr)
library(egg)

# Load Data
winSize = 5000

proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
data.dir = "/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP"

genomeWin = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/mm10_5kb-window.hSNPcounts",
                       header = F)
cnv = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/hdHetSNP/literature-Cahan-Alyriyami-Graubert_CNV.5kbWindow.hSNPcounts",
                 header = F)

# scatter theme
source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

# norm cnv to per 5kb
cnv = cnv %>% mutate(width = V3 - V2) %>% mutate(adjust = width / winSize) %>% mutate(norm_counts = V4 / adjust, tform = asinh(norm_counts))

# keep only standard chr in genomeWndows
genomeWin = genomeWin %>% filter(!grepl("random|chrUn|chrY|chrM", V1)) %>% mutate(tform = asinh(V4))
unique(genomeWin$V1)

# run KS
ks.res = ks.test(x = cnv$norm_counts, y = genomeWin$V4, exact = T)#, alternative = "less")
ks.res

# plot ecdf

p.ecdf = ggplot(cnv, aes(norm_counts)) +
  stat_ecdf(geom = "step", linewidth = 0.2, color = "red") +
  theme_classic() +
  xlab("HetSNPs per 5kb") +
  ylab("Cumulative distribution") +
  stat_ecdf(data = genomeWin, aes(V4), geom = "step",
            linewidth = 0.2, color = "black") + my.theme + xaxis.gap.fix #+ yaxis.gap.fix
p.ecdf

p.ecdf.egg = set_panel_size(p.ecdf, width = unit(1, "in"), height = unit(2,"in"))

ggsave(paste0(data.dir,"/",Sys.Date(), "_hetSNPs_WinGenomeVsWinCNV_KS.pdf"), # check which CNV data is used
       plot = p.ecdf.egg,
       device = "pdf",
       unit = "in",
       width = 2,
       height = 3,
       dpi = 300)

# percent hetSNP in peak
sum(cnv$V4) / sum(genomeWin$V4) * 100

ggplot(cnv, aes(asinh(V4))) +
  stat_ecdf(geom = "step", linewidth = 0.2, color = "red") +
  theme_classic() +
  xlab("HetSNPs per 5kb") +
  ylab("Cumulative distribution") +
  stat_ecdf(data = genomeWin, aes(asinh(V4)), geom = "step",
            linewidth = 0.2, color = "black") + my.theme + xaxis.gap.fix

