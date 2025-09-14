#bar plot for num NOD SNPs each chr
library(egg)
library(ggplot2)
library(dplyr)
library(stringr)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

snps_type = "homo"

if(snps_type == "all"){
  search_pattern = "PASS.hqHetSNP"
  plot_name = "All"
  axis_title = "All"
}else if(snps_type == "homo"){
  search_pattern = "onlyPASS"
  plot_name = "Homo"
  axis_title = "Homozygous"
}else if(snps_type == "het"){
  search_pattern = "chrUCSC.hqHetSNP.vcf"
  plot_name = "Het"
  axis_title = "Heterozygous"
}

snp_bed = list.files("/Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/vcf/by_chr",
                     full.names = T,
                     pattern = paste0("NOD_.*",search_pattern,".*bed"))
names(snp_bed) = basename(snp_bed) %>% str_match(pattern = "(chr[0-9X]*).bed") %>% .[,2]
snp_chr = lapply(snp_bed, FUN = function(x){rtracklayer::import(x) %>% length()}) %>% unlist()
snp_chr = data.frame(chromosome = names(snp_chr), num_snps = snp_chr)
#snp_chr = snp_chr[paste0("chr", c(seq_len(19),"X")),]
snp_chr$chromosome = factor(snp_chr$chromosome, levels = paste0("chr", c(seq_len(19),"X")))

ggplot(snp_chr, aes(x = chromosome, y = num_snps)) +
  geom_bar(stat = "identity",position = "dodge",width=0.75) +
  theme_classic() +
  my.theme + xaxis.gap.fix +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

chr.sizes = read.delim("/Scottbrowne/members/smd/genomes/mm10/mm10.chrom.sizes", header = F)
chr.sizes = chr.sizes[!grepl("random|chrUn|chrM|chrY", chr.sizes$V1),]
colnames(chr.sizes) = c("chromosome", "size")        

snp_chr = left_join(snp_chr, chr.sizes, by = "chromosome")
snp_chr$rate = snp_chr$num_snps / (snp_chr$size / 1e6)
snp_chr$chromosome = factor(snp_chr$chromosome, levels = paste0("chr", c(seq_len(19),"X")))


coeff = 100 # axis scaling

p.snp = ggplot(snp_chr, aes(x = chromosome, y = num_snps)) +
  geom_bar(stat = "identity",width=0.75, color = "grey30", fill = "grey30") +
  geom_line(data = snp_chr, aes(y = rate*coeff, x = chromosome), group = 1, stat = "identity", color = "red") +
 # geom_hline(yintercept = 1*coeff, linetype = "dashed", color = "orange") + 
  theme_classic() +
  my.theme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  scale_y_continuous(paste0("# ",axis_title," SNPs"), sec.axis = sec_axis( trans= ~ . / coeff, name=paste0(axis_title, " SNPs per Mb")), 
                     expand = expansion(mult = c(0, .05))) + 
  theme(
    axis.title.y = element_text(color = "grey30"),
    axis.title.y.right = element_text(color = "red")
  ) 

p.snp.egg = set_panel_size(p.snp, width = unit(3, "in"), height = unit(2, "in"))

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_",plot_name,"-snps-per-chr.pdf"),
       p.snp.egg,
       width = 4,
       height = 3,
       units = "in")
