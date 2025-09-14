library(ggplot2)
library(stringr)
library(dplyr)
library(egg)
library(GRaNIE)

genes =  read.delim("/Scottbrowne/seq/tmp/devoes/SD030/genes4QC.txt",header = F)
genes = paste0("^", genes$V1, "$")
#genes = c("Il11ra1","Il11ra2","Ccl19","Ccl27a","Ccl27b","Ccl21a","Ccl21b","Ccl21d")

proj_path = "/Scottbrowne/members/smd/Projects/SD030"
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                     "2024-10-06","_peak-gene_pearson_GRN.rds"))
group = GRN@data$metadata$archr_peak_group

labels = read.delim("/Scottbrowne/members/smd/Projects/SD030/sample_inputs/Group2Label.txt",
                    header = F)
labels = setNames(labels$V2, labels$V1)
sort_key = as.vector(labels)

group_label = labels[group]

rna = getCounts(GRN, type = "rna")
gene_transl = GRN@annotation$genes[,c("gene.ENSEMBL","gene.name")]
colnames(gene_transl) = c("ENSEMBL","gene.name")
rna = left_join(rna, gene_transl, by = "ENSEMBL")
rownames(rna) = rna$ENSEMBL
#rownames(rna) = rna$gene.name

rna = rna[grepl(paste0(genes,collapse = "|"),rna$gene.name) ,
          !(grepl("ENSEMBL",colnames(rna)))] 
rownames(rna) = rna$gene.name
rna = rna[,!(grepl("gene.name",colnames(rna)))]
rna = rna %>% t() %>% as.data.frame()
rna$group_label = group_label

#%>% t() %>% as_tibble()

rna_mean = rna %>% group_by(group_label) %>% summarise(across(.cols = everything(), mean))


rna_mean_mtx = as.matrix(rna_mean[,!grepl("group",colnames(rna_mean))])
rownames(rna_mean_mtx) = rna_mean$group_label
rna_mean_mtx = rna_mean_mtx[sort_key,]
rna_mean_mtx_center = scale(rna_mean_mtx, center = T, scale = F)

# heatmap ####

hm.genes = ComplexHeatmap::Heatmap(rna_mean_mtx_center,
                        cluster_rows = F,
                        show_row_names = T,
                        row_names_side = "left",
                        cluster_columns = T,
                       # column_order = match(colnames(rna_mean_mtx_center), genes %>% str_remove_all(pattern = "\\^|\\$")),
                        border = "black",
                        height = unit(4, "in"), 
                        width = unit(4, "in") 
                        )


source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")


pdf(paste0(proj_path,"/figures/",Sys.Date(),"QC_rna-expression.pdf"),
       width = 12,
       height =7)
hm.genes
dev.off()
