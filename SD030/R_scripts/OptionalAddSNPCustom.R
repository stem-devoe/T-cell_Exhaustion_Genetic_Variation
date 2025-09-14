# Add SNP data ####

# do custom add
snp_gr = lapply(c(seq_len(19),"X"),
                FUN = function(x, peak_gr){
                  snp_gr = rtracklayer::import(paste0("/Scottbrowne/members",
                                                      "/smd/genomes/NOD_SHILTJ_v5_REL1505",
                                                      "/vcf/by_chr",
                                                      "/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.PASS.hqHetSNP.vcf.chr",
                                                      x, ".bed"))
                  overlaps = findOverlaps(snp_gr, peak_gr) 
                  snp_gr[queryHits(overlaps)]
                },
                peak_gr) %>% 
  purrr::reduce(.,c) 
snp_gr$annotation = snp_gr$name

addSNPData_CUSTOM <- function(GRN, SNPs.gr){
  peaks.gr  = GRaNIE:::.constructGRanges(GRN@data$peaks$counts_metadata, 
                                         seqlengths = GRaNIE:::.getChrLengths(genomeAssembly), genomeAssembly)
  overlapsAll = GenomicRanges::findOverlaps(peaks.gr, SNPs.gr, 
                                            minoverlap = 1,
                                            type = "any", select = "all",
                                            ignore.strand = TRUE)
  
  
  
  query_row_ids   = S4Vectors::queryHits(overlapsAll)
  subject_row_ids = S4Vectors::subjectHits(overlapsAll)
  
  query_overlap_df     = as.data.frame(S4Vectors::elementMetadata(peaks.gr)[query_row_ids, "peakID"])
  subject_overlap_df   = as.data.frame( SNPs.gr)[subject_row_ids, c("seqnames", "start", "annotation")]
  
  overlaps.df = cbind.data.frame(query_overlap_df,subject_overlap_df) %>%
    dplyr::mutate(seqnames = as.character(.data$seqnames)) %>% tibble::as_tibble()
  colnames(overlaps.df) = c("peakID", "SNP_chr", "SNP_start", "SNP_rsid")
  
  
  # Summarize on the peak level before integrating
  overlaps.grouped.df = overlaps.df %>%
    dplyr::group_by(.data$peakID) %>%
    dplyr::summarise(SNP_n = dplyr::n(), 
                     SNP_rsid = paste0(.data$SNP_rsid, collapse = ","), 
                     SNP_start = paste0(.data$SNP_start, collapse = ",")) #, SNP_origin = paste0(.data$association, collapse = ",")) 
  
  GRN@data$peaks$counts_metadata = dplyr::left_join(GRN@data$peaks$counts_metadata, overlaps.grouped.df, 
                                                    by = "peakID", multiple = "all")
  GRN@data$peaks$counts_metadata$SNP_n[which(is.na(GRN@data$peaks$counts_metadata$SNP_n))] = 0
  
  futile.logger::flog.info(paste0(" Added SNP information to GRN@data$peaks$counts_metadata"))
  
  futile.logger::flog.info(paste0(" Statistics: "))
  for (i in sort(unique(GRN@data$peaks$counts_metadata$SNP_n))) {
    futile.logger::flog.info(paste0("  Number of peaks with ", i, " SNPs associated: ", dplyr::filter(GRN@data$peaks$counts_metadata, .data$SNP_n == i) %>% nrow()))
  }
  
  
  # GRaNIE:::.printExecutionTime(start, prefix = "")
  GRN
  
}

GRN = addSNPData_CUSTOM(GRN, snp_gr)

saveRDS(GRN, paste0(proj_path,"/data/GRANIE/",
                    Sys.Date(),"_addedSNPData_GRN.rds"))

#GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
#                     "2024-07-11","_addedSNPData_GRN.rds")) 2024-07-10 may have mm39 in gene annos

