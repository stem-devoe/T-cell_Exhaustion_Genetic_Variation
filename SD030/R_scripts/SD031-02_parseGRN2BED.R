save_date = "2024-10-06"
proj_path = "/Scottbrowne/members/smd/Projects/SD030"
library(dplyr)
library(stringr)

GRN_connections.filt = readRDS(paste0(proj_path, "/data/GRANIE/",
                                      save_date, "_GRN-connections_TFpeak0.3FDR_peakgene0.2FDR_pearson.rds"))



p2g = readRDS(paste0(proj_path,"/data/GRANIE/",
                     save_date,"_peak2gene_pearson.rds"))
# when does granie do mult testing corrections
p2g[p2g$gene.ENSEMBL == "ENSMUSG00000027087",c(1,9)]
GRN_connections.filt[GRN_connections.filt$gene.ENSEMBL == "ENSMUSG00000027087", c(4,17,18)]

# add mult testing corrections
p2g$peak_gene.p_adj = p.adjust(p2g$peak_gene.p_raw, method = "fdr")

# add gene name to p2g
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                     save_date,"_TF-peak-gene_pearson_GRN.rds"))

p2g = left_join(p2g, GRN@annotation$genes[,c('gene.ENSEMBL','gene.name')], by = "gene.ENSEMBL")

# format p2g for bedtools intersect ####
p2g.tbl = with(p2g, tibble(peak = peak.ID,
                           target = gene.name,
                           distance = peak_gene.distance,
                           peak_gene.r = peak_gene.r,
                           peak_gene.p_raw = peak_gene.p_raw,
                           peak_gene.p_adj = peak_gene.p_adj))
p2g.tbl = p2g.tbl %>% mutate(chr = str_remove(peak, ":.*"), 
                             start = str_remove(str_remove(peak, ".*:"),"-.*"), # granie ID is 0-based
                             stop = str_remove(peak, ".*-")) %>%
  mutate(ID = paste0(chr,"_", start)) %>%
  relocate(chr) %>% relocate(start, .after = chr) %>% relocate(stop, .after = start) %>% relocate(ID, .after = stop)

write.table(p2g.tbl, paste0(proj_path, "/data/GRANIE/",
                            save_date,
                            "_peak2gene_pearson.bed"),
            row.names = F,col.names = F,quote = F, sep = "\t")


# format GRN_connections for bedtools intersect ####

# want
# chr start stop ID TF targetGene TF_peak.r TFpeak.fdr peak_gene.r peak_gene.fdr peak_gene.distance TF_gene.r
grn.tbl = with(GRN_connections.filt, tibble(peak = peak.ID, # granie's are 0-based
                    motif = TF.ID,
                    target = gene.name,
                    TF_peak.r = TF_peak.r,
                    TF_peak.fdr = TF_peak.fdr,
                    peak_gene.r = peak_gene.r,
                    peak_gene.p_adj = peak_gene.p_adj,
                    distance = peak_gene.distance,
                    TF_gene.r = TF_gene.r))

grn.tbl = grn.tbl %>% mutate(chr = str_remove(peak, ":.*"), 
                   start = str_remove(str_remove(peak, ".*:"),"-.*"), # granie ID is 0-based
                   stop = str_remove(peak, ".*-")) %>%
  mutate(ID = paste0(chr,"_", start)) %>%
  relocate(chr) %>% relocate(start, .after = chr) %>% relocate(stop, .after = start) %>% relocate(ID, .after = stop)

write.table(grn.tbl, paste0(proj_path, "/data/GRANIE/",
                            save_date,
                            "_GRN-connections_TFpeak0.3FDR_peakgene0.2FDR_pearson.bed"),
            row.names = F,col.names = F,quote = F, sep = "\t")
