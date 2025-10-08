#################################################################################################################
###########          Do NOT use this approach                                              ######################
###########          biomaRt will grab multiple ENSMBL Gene Ids for some TF                ######################
###########     create the table manually from HOCOMOCO database annotation file and MGI   ######################
###########        https://www.informatics.jax.org/batch/summary                           ######################
#################################################################################################################



# mgi to enzembl for hocomoco
# make sure you pull mm10
hoco_meta = "/Scottbrowne/members/smd/motifs/HOCOMOCO/v11/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv"
library(biomaRt)
hoco_meta = read.delim(hoco_meta)
#listEnsemblArchives()
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", 
                      version = 102, host = "https://nov2020.archive.ensembl.org") #, host =)
#listAttributes()
#listFilters(ensembl)# %>% grep(pattern = "mgi|ensembl", ignore.case = T, value = T)
#searchFilters(mart = ensembl, pattern = "ensembl.*id")
#searchAttributes()
motif_ids = getBM(attributes = c('mgi_id', 'ensembl_gene_id'),
      filters = "mgi_id",
      values = paste0("MGI:",hoco_meta$MGI), # need to format for use in filters MGI:<ID#>
      mart = ensembl)
trans_table = hoco_meta[,c('Transcription.factor','MGI','Model')]
trans_table$mgi_id = paste0("MGI:",trans_table$MGI)
trans_table = left_join(trans_table,motif_ids, by = "mgi_id")
trans_table = trans_table[,c('Transcription.factor','ensembl_gene_id','Model')]
colnames(trans_table) = c('SYMBOL','ENSEMBL','HOCOID')
head(trans_table)
write.table(trans_table,
          "/Scottbrowne/members/smd/motifs/HOCOMOCO/v11/GRaNIE/translationTable.csv",
          row.names = F,
          quote = F,
          sep = "\t")
