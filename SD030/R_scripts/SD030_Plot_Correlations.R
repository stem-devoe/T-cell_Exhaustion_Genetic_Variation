library(ggplot2)
library(stringr)
library(dplyr)
library(egg)
library(GRaNIE)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

fill.colors = c(
  Terminal = "#920000",
  Progenitor = "#009292",
  PD1loTIM3lo_TILs = "#000000",
  Memory = "#490092",
  Effector = "#db6d00",
  Naive = "#b66dff"
)

proj_path = "/Scottbrowne/members/smd/Projects/SD030"
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                     "2024-10-06","_peak-gene_pearson_GRN.rds"))
group = GRN@data$metadata$archr_peak_group

labels = read.delim("/Scottbrowne/members/smd/Projects/SD030/sample_inputs/Group2Label.txt",
                    header = F)
labels = setNames(labels$V2, labels$V1)
sort_key = as.vector(labels)

group_label = labels[group]


p2g.pairs = read.delim("/Scottbrowne/seq/tmp/devoes/SD030/p2g_peak_gene.txt")

# ATAC ####

peaks = p2g.pairs$peakID

atac = getCounts(GRN, type = "peaks")
rownames(atac) = atac$peakID %>% str_remove(pattern = "-.*") %>% 
  str_replace(pattern = ":", replacement = "_")

atac = atac[peaks ,
            !(grepl("peakID",colnames(atac)))] 
atac = atac %>% t() %>% as.data.frame()
#atac$group_label = group_label
atac$sample = rownames(atac)

# RNA ####

genes = paste0("^",p2g.pairs$gene,"$")

rna = getCounts(GRN, type = "rna")
gene_transl = GRN@annotation$genes[,c("gene.ENSEMBL","gene.name")]
colnames(gene_transl) = c("ENSEMBL","gene.name")
rna = left_join(rna, gene_transl, by = "ENSEMBL")
rownames(rna) = rna$ENSEMBL

rna = rna[grepl(paste0(genes,collapse = "|"),rna$gene.name) ,
          !(grepl("ENSEMBL",colnames(rna)))] 
rownames(rna) = rna$gene.name
rna = rna[,!(grepl("gene.name",colnames(rna)))]
rna = rna %>% t() %>% as.data.frame()
rna$group_label = group_label
rna$Phenotype = GRN@data$metadata$Phenotype
rna[c("SB2019_ENDO_Naive_1","SB2019_ENDO_Naive_3"),"Phenotype"] = "PD1loTIM3lo_TILs" # replace non exh TILs Pheno because I had it incorrect
rna$Study = GRN@data$metadata$Publication
rna$sample = rownames(rna)

# format ####

df = left_join(rna, atac, by = "sample")

# plot correlation ####
for(i in seq_len(nrow(p2g.pairs))){
  cor.test(df[,p2g.pairs$peakID[i]], df[,p2g.pairs$gene[i]], method = "pearson")
  p = df %>% rename(peak = p2g.pairs$peakID[i], gene = p2g.pairs$gene[i]) %>% #head(n = 1)
    ggplot( mapping = aes(x = peak, y = gene)) +
    geom_point(mapping = aes(fill = Phenotype, color = Phenotype, shape = Study), size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = F, linetype = 'dashed', color = "grey50") + 
    theme_classic() + my.theme + 
    xlab(paste0("Relative Accessibility\n",p2g.pairs$peakID[i])) + # corrected 20241031 (labels were flipped)
    ylab(paste0("Relative Expression\n", p2g.pairs$gene[i])) + # corrected 20241031
    scale_fill_manual(values = fill.colors) + scale_color_manual(values = fill.colors) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) # expand ranges so dots aren't sitting on edge of plot 
  p
  p.egg = set_panel_size(p, width = unit(2,"in"), height = unit(2, "in"))
  ggsave(paste0(proj_path,"/figures/",Sys.Date(),"_",
                p2g.pairs$gene[i], "_", p2g.pairs$peakID[i], 
                "_correlation_scatter.pdf"),
         plot=p.egg,
         width = 5,
         height = 5,
         device = "pdf")
}

