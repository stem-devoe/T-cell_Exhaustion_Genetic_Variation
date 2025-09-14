peaks = read.delim(paste0("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda",
                          "/150bp/All_Samples.fwp.filter.non_overlapping.renamed.bed"),
                   header = F)
# seqname, source, feature, start, end, score, strand, frame, attribute
gff = data.frame(seqname = peaks[,1],
                 source = "SD029",
                 feature="RegulatoryElement",
                 start = peaks[,2] + 1, # gff is 1-based
                 end = peaks[,3],
                 score = peaks[,5],
                 strand = ".",
                 frame = 0,
                 attribute = paste0("ID=",peaks[,4])
)
write.table(gff,
            col.names = F,
            row.names = F,
            sep = "\t",
            file = paste0("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda",
                          "/150bp/All_Samples.fwp.filter.non_overlapping.renamed.gff"),
            quote = F)
