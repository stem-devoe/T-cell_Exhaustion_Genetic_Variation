save_date = "2024-10-06"
proj_path = "/Scottbrowne/members/smd/Projects/SD030"
library(dplyr)
library(stringr)


p2g = readRDS(paste0(proj_path,"/data/GRANIE/",
                     save_date,"_peak2gene_pearson.rds"))

# add mult testing corrections
p2g$peak_gene.p_adj = p.adjust(p2g$peak_gene.p_raw, method = "fdr")

# add gene name to p2g
GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
                     save_date,"_TF-peak-gene_pearson_GRN.rds"))

p2g = left_join(p2g, GRN@annotation$genes[,c('gene.ENSEMBL','gene.name')], by = "gene.ENSEMBL")

mart = biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "102")
tss <- biomaRt::getBM(attributes = c("transcription_start_site", 
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id", "ensembl_transcript_id",
                            "ccds", 'transcript_gencode_basic','external_transcript_name',
                            "ensembl_gene_id", "external_gene_name"),
             filters = "ensembl_gene_id", values = unique(p2g$gene.ENSEMBL),
             mart = mart)
tss.filt = tss[tss$ccds != "",] 
tss.filt = tss.filt[order(tss.filt$external_transcript_name), ] # sort by transcript name in ensembl (NOT ID); I'm fairly sure the GENE-# with the lowest number is the canonical one
tss.filt.uni = tss.filt[!duplicated(tss.filt$ensembl_gene_id),] # take first to get lowest num in transscript NAME
  
p2g$ensembl_gene_id = p2g$gene.ENSEMBL
p2g = left_join(p2g, tss.filt.uni, by = "ensembl_gene_id")


# format p2g for UCSC interact track ####

p2g.tbl = with(p2g, tibble(chrom = str_remove(peak.ID, ":.*"),
                           chromStart = 0,
                           chromEnd = 0,
                           name = paste0(peak.ID,"_",gene.name),
                           score = 0,
                           value = peak_gene.r,
                           exp = "SD030_GRaNIE",
                           color = 0,
                           sourceChrom = str_remove(peak.ID, ":.*"),
                           sourceStart = str_remove(str_remove(peak.ID, ".*:"),"-.*"), # granie ID is 0-based
                           sourceEnd = str_remove(peak.ID, ".*-"),
                           sourceName = peak.ID,
                           sourceStrand = ".",
                           targetChrom = gene.chr,
                           targetStart = transcription_start_site,
                           targetEnd = transcription_start_site + 1,
                           targetName = gene.name,
                           targetStrand = gene.strand))
p2g.tbl$chromStart = p2g.tbl$sourceStart
p2g.tbl$chromEnd = p2g.tbl$sourceEnd

chromStart_bool = which(p2g.tbl$targetStart < p2g.tbl$sourceStart)
chromEnd_bool = which(p2g.tbl$targetEnd > p2g.tbl$sourceEnd)

p2g.tbl$chromStart[chromStart_bool] = p2g.tbl$targetStart[chromStart_bool]
p2g.tbl$chromEnd[chromEnd_bool] = p2g.tbl$targetEnd[chromEnd_bool]

# add color

fill.colors = hcl.colors(201,palette = "RdBu", rev = F)
p2g.tbl$color = fill.colors[round(p2g.tbl$value * 100)+101]

# save

# write.table(p2g.tbl, paste0(proj_path, "/data/GRANIE/",
#                             save_date,
#                             "_peak2gene_pearson_interact_UCSC.txt"),
#             row.names = F,col.names = F,quote = F, sep = "\t")
# 
# # make scale
# v.scale = data.frame(x = seq_len(201), y = seq_len(201), fill = fill.colors)
# p.scale = ggplot(v.scale, aes(y=x,x = y,fill=x)) + geom_point() + scale_fill_distiller(palette =  "RdBu",direction = 1 )
# ggsave("/Scottbrowne/members/smd/Projects/SD030/figures/loops_scale.pdf",
#        plot = p.scale,
#        device = "pdf",
#        width = 4,
#        height = 4)

# save selected genes

genes = c('Cxcr4',
          'Apoa2',
          'Slamf6',
          'Usp44',
          'Tbx21',
          'Cd47',
          'Dusp5',
          'Itgav',
          'Ptpn22',
          'Spry1',
          'Il21',
          'Id3',
          'Bach2',
          'Ptpn11',
          'Lag3',
          'Mt3',
          'Mt1',
          'Mt2',
          'Eomes')
p2g.tbl.select = p2g.tbl %>% filter(targetName %in% genes)
write.table(p2g.tbl.select, paste0(proj_path, "/data/GRANIE/",
                           save_date,
                           "_SELECT_peak2gene_pearson_interact_UCSC.txt"),
           row.names = F,col.names = F,quote = F, sep = "\t")
