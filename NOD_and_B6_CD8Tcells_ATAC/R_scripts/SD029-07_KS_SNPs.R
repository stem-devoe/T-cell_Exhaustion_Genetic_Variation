#stat_ecdf
#ecdf

library(ggplot2)
library(dplyr)
library(stringr)
library(egg)

# Load Data

proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
fold_change = 1.5
half_width = 150
padj = 0.1
pheno = c("CART","LCMVcl13","Naive")

# scatter theme
source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

# color palettes 
sort.key = c("B6_Naive",
             "NOD_Naive",
             "B6_Memory",
             "NOD_Memory",
             "B6_LCMVcl13",
             "NOD_LCMVcl13",
             "B6_CART",
             "NOD_CART")

fill.colors = c(B6_Naive = "#ff33ff",
                NOD_Naive = "#cc99ff",
                B6_Memory = "#ff6600",
                NOD_Memory = "#ff9999",
                B6_LCMVcl13 = "#009900",
                NOD_LCMVcl13 = "#99cc66",
                #B6_CART = "#3333ff",
                B6_CART = "#4969ff",
                NOD_CART = "#6699ff")

# sig peaks per comparison
peak.compare = list.files(paste0(proj.dir, "/data/",half_width, "bp/padj_",padj,"/sig_peaks/FC_",fold_change),
                          pattern = "up.bed|down.bed",
                          full.names = T) %>% 
  grep(pattern = "NOD-B6", value = T) %>%
  grep(pattern = paste0(pheno, collapse = "|"), value = T)


# all peaks with SNPs per peak
all.snp.counts = read.delim(paste0(proj.dir, "/data/All_Samples.fwp.filter.non_overlapping.renamed.countPASSandHetSNPs.bed"), header = F)
colnames(all.snp.counts) = c("chr","start","stop","ID","score","strand","SNPs")

homo.snp.counts = read.delim(paste0(proj.dir, "/data/All_Samples.fwp.filter.non_overlapping.renamed.countOnlyPASS.bed"), header = F)
colnames(homo.snp.counts) = c("chr","start","stop","ID","score","strand","SNPs")

het.snp.counts = left_join(all.snp.counts, homo.snp.counts[,c('ID','SNPs')], by = "ID")
colnames(het.snp.counts) = c("chr","start","stop","ID","score","strand","all.SNPs","homo.SNPs")
het.snp.counts$het.SNPs = with(het.snp.counts, all.SNPs - homo.SNPs)
quantile(het.snp.counts$het.SNPs, probs = seq(0.985,1,0.001))
length(which(het.snp.counts$het.SNPs > 0 )) # 947 peaks with hSNPs of ~69000. sub 1%


# 0. define background peaks ####

dt.res = readRDS(paste0(proj.dir, "/data/",half_width, "bp/padj_",padj,"/decideTests_FC",fold_change,"_countsFiltered.rds"))

all.sig.peaks = lapply(peak.compare, FUN = function(x){x = read.delim(x, header = F); x$V4}) %>% 
  unlist() %>% unique()
bg.peaks = data.frame(ID = rownames(dt.res)[!(rownames(dt.res) %in% all.sig.peaks)])
bg.snps = left_join(bg.peaks, het.snp.counts, by = "ID") %>% 
  relocate(ID, .after = stop)


#  pairwise ####

# get SNPs df for significant peaks
up.snps = lapply(pheno, FUN = getSNPs, 
                 direction = "up",
                 sigPeakFiles = peak.compare, 
                 snpCounts = het.snp.counts)
names(up.snps) = paste0("NOD ", pheno, " > B6 ", pheno)
down.snps = lapply(pheno, FUN = getSNPs, direction = "down", 
                   sigPeakFiles = peak.compare, 
                   snpCounts = het.snp.counts)
names(down.snps) = paste0("NOD ", pheno, " < B6 ", pheno)

sig.snps = c(up.snps,down.snps)

# peaks with SNP 
nrow(bg.snps)
table(bg.snps$homo.SNPs > 0)
table(bg.snps$all.SNPs > 0)
table(bg.snps$het.SNPs > 0)
names(sig.snps)
lapply(sig.snps, nrow)
lapply(sig.snps, FUN = function(x){table(x$homo.SNPs > 0)})
lapply(sig.snps, FUN = function(x){table(x$all.SNPs > 0)})
lapply(sig.snps, FUN = function(x){table(x$het.SNPs > 0)})

# average num SNPs per set
bg.snps$homo.SNPs %>% mean
bg.snps$all.SNPs %>% mean
bg.snps$het.SNPs %>% mean
lapply(sig.snps, function(x){mean(x$homo.SNPs)})
lapply(sig.snps, function(x){mean(x$all.SNPs)})
lapply(sig.snps, function(x){mean(x$het.SNPs)})

# run KS test
all.ks.res = lapply(sig.snps, FUN = runKS,  bgSNPs = bg.snps, snps_type = "all")
homo.ks.res = lapply(sig.snps, FUN = runKS,  bgSNPs = bg.snps, snps_type = "homo")
het.ks.res = lapply(sig.snps, FUN = runKS,  bgSNPs = bg.snps, snps_type = "het")
# be sure check pvalue

# plot eCDF ####

p.all =  plotECDF(sig.snps,
                      fillColors = fill.colors, 
                      bgSNPs = bg.snps, 
                      snps_type = "all",
                  linewidth = 0.35,
                      customTheme = my.theme)
p.all = p.all + xaxis.gap.fix + yaxis.gap.fix
p.all
p.all.egg = set_panel_size(p.all, width = unit(1, "in"), height = unit(2,"in"))

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_allSNPS_ecdf.pdf"),
       plot = p.all.egg,
       device = "pdf",
       unit = "in",
       width = 2,
       height = 3,
       dpi = 300)

p.homo =  plotECDF(sig.snps,
                  fillColors = fill.colors, 
                  bgSNPs = bg.snps, 
                  snps_type = "homo",
                  linewidth = 0.35,
                  customTheme = my.theme)
p.homo = p.homo + xaxis.gap.fix + yaxis.gap.fix
p.homo
p.homo.egg = set_panel_size(p.homo, width = unit(1, "in"), height = unit(2,"in"))
ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_homoSNPS_ecdf.pdf"),
       plot = p.homo.egg,
       device = "pdf",
       unit = "in",
       width = 2,
       height = 3,
       dpi = 300)

p.het =  plotECDF(sig.snps,
                  fillColors = fill.colors, 
                  bgSNPs = bg.snps, 
                  snps_type = "het",
                  linewidth = 0.35,
                  customTheme = my.theme)
p.het = p.het + xaxis.gap.fix + yaxis.gap.fix
p.het
p.het.egg = set_panel_size(p.het, width = unit(1, "in"), height = unit(2,"in"))
ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_hetSNPS_ecdf.pdf"),
       plot = p.het.egg,
       device = "pdf",
       unit = "in",
       width = 2,
       height = 3,
       dpi = 300)


## b. bar plot of sig peak counts #####

dt.res.summary = apply(dt.res, MARGIN = 2, FUN = table) %>% as.data.frame()
rownames(dt.res.summary) = c("Down","NotSig","Up")
dt.res.summary = dt.res.summary[,grepl(paste0(pheno,collapse = "|"), colnames(dt.res.summary))]
dt.res.summary$Change = rownames(dt.res.summary)
dt.res.summary = reshape2::melt(dt.res.summary)
dt.res.summary$variable = str_remove(dt.res.summary$variable, "NOD-B6_")
dt.res.summary$Change = str_replace(dt.res.summary$Change, pattern = "Down", replacement = "NOD < B6") %>%
  str_replace(pattern = "Up", replacement = "NOD > B6")
dt.res.summary = dt.res.summary[!grepl("NotSig",dt.res.summary$Change),]
dt.res.summary
p.bar = ggplot(dt.res.summary, aes(x = variable, y = value, fill = Change))  + 
  geom_bar(position="stack", stat="identity", width = 0.75) + 
  theme_classic() + 
  my.theme + xaxis.gap.fix +
  scale_fill_manual(values=c(`NOD < B6` = "grey20", `NOD > B6` = "grey60")) + ylab("number of peaks") +
  theme(legend.position = "none") + xlab("")
p.bar

p.bar.egg = set_panel_size(p.bar, width = unit(1, "in"), height = unit(2,"in"))

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(), "_bar.pdf"),
       plot = p.bar.egg,
       device = "pdf",
       unit = "in",
       width = 2,
       height = 3,
       dpi = 300)


# Functions #####

plotECDF <- function(sigSNPs, numCompare, fillColors, bgSNPs, snps_type,
                     customTheme, linewidth = 0.2){
  bgSNPs = data.frame(SNPs = bgSNPs[,paste0(snps_type, ".SNPs")])
  p.all = ggplot(bgSNPs, aes(SNPs)) + 
    stat_ecdf(geom = "step", linewidth = linewidth) + 
    theme_classic( ) + 
    customTheme +
    xlab("SNPs/peak") + 
    ylab("Cumulative distribution")
  for(i in seq_along(sigSNPs)){
    comparison = names(sigSNPs)[i]
    phenotype = str_match(comparison, pattern = "NOD (.*) (<|>)") %>% .[,2]
    if(grepl("<",comparison)){
      strain = "B6"
    } else{
      strain = "NOD"
    }
    sig.df = data.frame(SNPs = sigSNPs[[i]][,paste0(snps_type,".SNPs")])
    p.all = p.all + stat_ecdf(data = sig.df,
                                      geom = "step", 
                                      linewidth = linewidth,
                                      color = fill.colors[paste0(strain,"_",phenotype)])
  }
  return(p.all)
}

runKS <- function(sigSNPs, bgSNPs, snps_type, runAsExact = T){
  sigSNPs = sigSNPs[,paste0(snps_type,".SNPs")] %>% as.vector()
  bgSNPs = bgSNPs[,paste0(snps_type,".SNPs")] %>% as.vector()
  ks.res = ks.test(x = sigSNPs, 
                   y = bgSNPs,
                   exact = runAsExact)
  return(ks.res)
}

getSNPs <- function(phenotype, direction = c("up","down"), sigPeakFiles, snpCounts){
  peak.file = grep(pattern = phenotype, sigPeakFiles, value = T) %>% 
    grep(pattern = direction, x = ., value = T)
  peak = read.delim(peak.file, header = F)
  colnames(peak) = c("chr","start","stop","ID")
  snps = left_join(peak, snpCounts, by = c("chr","start","stop", "ID"))
  #snps = snpCounts[snpCounts$ID %in% peak[,4], ]
  return(snps)
}
