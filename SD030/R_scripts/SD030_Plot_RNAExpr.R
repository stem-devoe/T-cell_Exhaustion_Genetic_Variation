library(ggplot2)
library(stringr)
library(dplyr)
library(egg)
library(GRaNIE)

genes = c(
  "Pdcd1","Ifng","Sell","Tcf7",#"Havcr2","Gzmb","Tnf","Prf1",
  "Eomes","Tnfrsf9","Id3","Il21","Cxcr4","Ptpn7"
  #"Itgav",
   #       "Havcr2",
    #      "Itk",
          #"Lck",
         # "Il27ra",
     #     "Cxcr4",
          #'Pdcd1','Tcf7','Ifng','Tnf', # sanity checks
         # "Ptpn7",
         # "Ptpn18",
         # "Ptprv",
      #   "Cd244a"
      #    "Cd48", 
        #  "Socs4"
          )
genes = paste0("^",genes,"$")

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
sort_order = lapply(sort_key, 
                    FUN = function(key,groups){grep(key,groups, 
                                                    ignore.case = T, 
                                                    value = T)},
                      rownames(rna_mean_mtx)) %>% unlist()
rna_mean_mtx = rna_mean_mtx[sort_order,]
rna_mean_mtx_center = scale(rna_mean_mtx, center = T, scale = F)


# dot bar graph ####

rna.dot.df = rna_mean_mtx_center %>% as.data.frame() %>% mutate(group = rownames(rna_mean_mtx_center))
rna.dot.df = reshape2::melt(rna.dot.df)
rna.dot.df$group = factor(rna.dot.df$group, levels = rev(sort_order))

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

p.bardot = ggplot(rna.dot.df, aes(x = value, y = group, group = variable)) + 
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey40") + 
  geom_point(size = 2, col = "black") + facet_grid(. ~ variable, scales = "free") +
  theme_classic() + my.theme + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) # expand ranges so dots aren't sitting on edge of polot 
p.bardot

p.egg = set_panel_size(p.bardot, width = unit(45, "mm"), height = unit(183.7, "mm"))

ggsave(paste0(proj_path,"/figures/rna-expression2.pdf"),
       plot = p.egg,
       width = 24,
       height = 12,
       device = "pdf")

