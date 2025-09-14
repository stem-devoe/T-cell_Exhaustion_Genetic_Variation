library(dplyr)
library(stringr)
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicInteractions)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)


# NEED TO DECIDE NOD SNPs TYPE
# currently only does homozygous

# sizes width = 4 height = 5 is good if only CART and Naive
# width = 4 and height = 5.5 for CART, LCMV, and Naive

# parameters ####

regions = read.delim("/Scottbrowne/members/smd/Projects/SD029/plotRegions/eGRN_regions.txt")


for(i in seq_len(nrow(regions))){
  genome = "mm10"
  chr = regions$chr_view[i] # "chr2"
  gene = regions$Gene[i] # "Itgav"
  pos.start = regions$start_view[i] #83693000  
  pos.end = regions$end_view[i] # 83750000 
  multiome.bw.lim = c(0,regions$bw_ref_max[i])
  gv.bw.lim = c(0,regions$bw_gv_max[i])
  highlight.start = regions$start_highlight[i]	# 83694022
  highlight.end = regions$end_highlight[i]
  motif = regions$motif_disrupt[i]
  regulon=regions$Regulon[i]
  loop1 = regions$Loop1[i]
  loop2 = regions$Loop2[i]
  
  
  fill.colors = c(B6_Naive = "#ff33ff",
                  NOD_Naive = "#cc99ff",
                  B6_Memory = "#ff6600",
                  NOD_Memory = "#ff9999",
                  B6_LCMVcl13 = "#009900",
                  NOD_LCMVcl13 = "#99cc66",
                  #B6_CART = "#3333ff",
                  B6_CART = "#4969ff",
                  NOD_CART = "#6699ff")
  
  # paths ####
  
  gv_meta = read.delim("/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv")
  
  # bigwigs
  gv_bw_files = list.files("/scratch2/devoes/SD029/bigwig_frip_groupAvg",
                           full.names = T,
                           pattern = ".bw")# %>% grep(pattern = "Naive|CART", value = T)
 # gv_bw_files = c(gv_bw_files[grepl("CART", gv_bw_files)], gv_bw_files[grepl("Naive", gv_bw_files)])
  gv_bw_files = c(gv_bw_files[grepl("CART", gv_bw_files)],
                  gv_bw_files[grepl("LCMV", gv_bw_files)],
                  gv_bw_files[grepl("Naive", gv_bw_files)])
  names(gv_bw_files) = gv_bw_files %>% basename() %>% str_remove(pattern = "_frip.*")
  
  multiome_bw_files = list.files("/Scottbrowne/members/smd/Projects/SD031/scenicplus/B16_CART/ATAC/consensus_peak_calling/pseudobulk_bw_files",
                                 full.names = T,
                                 pattern = ".bw")
  names(multiome_bw_files) = basename(multiome_bw_files) 
  multiome_bw_files = multiome_bw_files[c(paste0(seq_along(multiome_bw_files) - 1 , ".bw"))]
  
  # peaks
  gv_consensus = rtracklayer::import("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.renamed.bed")
  gv_consensus = gv_consensus[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                               ranges = IRanges(start = pos.start, end = pos.end),
                                                               strand = "*"),
                                                       gv_consensus)),]
  
  gv_cart_acquired = readRDS("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/150bp/padj_0.1/FC_1.5/CART/Gained_CART-to-Naive_any-strain_peak_subset.rds")
  gv_cart_acquired = gv_cart_acquired[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                                       ranges = IRanges(start = pos.start, end = pos.end),
                                                                       strand = "*"),
                                                               gv_cart_acquired)),]
  
  gv_cart_nodloss = rtracklayer::import("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5/NOD-B6_CART.down.bed")
  gv_cart_nodloss = gv_cart_nodloss[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                                     ranges = IRanges(start = pos.start, end = pos.end),
                                                                     strand = "*"),
                                                             gv_cart_nodloss)),]
  
  multiome_consensus = rtracklayer::import("/Scottbrowne/members/smd/Projects/SD031/scenicplus/B16_CART/ATAC/consensus_peak_calling/consensus_regions.bed")
  multiome_consensus = multiome_consensus[grepl(chr,seqnames(multiome_consensus))]
  seqlevels(multiome_consensus) = seqlevels(multiome_consensus) %>% grep(pattern = "chr", value = T)
  multiome_consensus = multiome_consensus[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                                           ranges = IRanges(start = pos.start, end = pos.end),
                                                                           strand = "*"),
                                                                   multiome_consensus)),]
  
  # polymorphisms
  nod_snp = rtracklayer::import(paste0("/Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/vcf/by_chr/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.onlyPASS.vcf.", chr, ".bed"),
                                extraCols = c())
  nod_snp = nod_snp[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                     ranges = IRanges(start = pos.start, end = pos.end),
                                                     strand = "*"),
                                             nod_snp)),]
  
  # motifs
  motif_gr = read.delim("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_consensus_b6sites.bed",
                        header = F)
  motif_gr = motif_gr[motif_gr$V8 %in% motif, c(1,2,3,8)]
  colnames(motif_gr) = c("chr","start","end","motif")
  motif_gr = makeGRangesFromDataFrame(motif_gr, starts.in.df.are.0based = T, keep.extra.columns = T)
  motif_gr = motif_gr[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                       ranges = IRanges(start = pos.start, end = pos.end),
                                                       strand = "*"),
                                               motif_gr)),]
  # multiome tracks ####
  
  multiome.bw.dtrack = lapply(multiome_bw_files,
                              FUN = DataTrack,
                              genome = "mm10",
                              type = "h", 
                              window = -1, # add smoothing effect
                              fill = "grey20",
                              chromosome = chr, 
                              name = "bigwig",
                              ylim = multiome.bw.lim,
                              yTicksAt= multiome.bw.lim[2],
                              fontsize.legend = 6,
                              col = "grey50",
                              col.histogram=NA)
  
  multiome.consensus.atrack = AnnotationTrack(multiome_consensus,
                                              genome = "mm10",
                                              chromosome = chr,
                                              start = pos.start,
                                              end = pos.end,
                                              col = "black",fill = "black",
                                              name = "consensus peaks", 
                                              collapse= T,
                                              shape = "box")
  
  peak_gr = GRanges(Rle(chr),IRanges(start = loop1, end = loop1),
                    Rle(strand(c("*"))))
  tss_gr = GRanges(Rle(chr),IRanges(start = loop2, end = loop2),
                   Rle(strand(c("*"))))
  loop_gr = GenomicInteractions(peak_gr,tss_gr,counts = 10)
  egrn.track = InteractionTrack(name = "eGRN_peak-to-gene", loop_gr, chromosome = chr)
  
  # NOD v B6 tracks ####
  
  # bigwigs
  
  gvDataTrack = function(bw_file, fill_colors, chr, bw.lim){
    DataTrack(bw_file, genome = "mm10",
              type = "h",
              window = -1, # add smoothing effect
              chromosome = chr,
              ylim = bw.lim,
              yTicksAt = bw.lim[2],
              fontsize.legend = 6,
              col = fill.colors[bw_file %>% basename() %>% str_remove(pattern = "_frip.*") %>% str_remove(pattern = "_Exh")])
  }
  gv.bw.dtrack = lapply(gv_bw_files, FUN = gvDataTrack, fill.colors, chr, gv.bw.lim)
  
  
  # peaks
  gv.consensus.atrack = AnnotationTrack(gv_consensus,
                                        genome = "mm10",
                                        chromosome = chr,
                                        start = pos.start,
                                        end = pos.end,
                                        col = "black", fill = "black",
                                        name = "consensus peaks", 
                                        collapse = T,
                                        shape = "box")
  # loss.atrack = AnnotationTrack(gv_cart_nodloss,
  #                               genome = "mm10",
  #                               chromosome = chr,
  #                               start = pos.start,
  #                               end = pos.end,
  #                               col = "purple",
  #                               name = "NOD < B6",
  #                               shape = "box")
  
  
  # motif tracks ####
  
  motif.atrack = AnnotationTrack(motif_gr,
                                 genome = "mm10",
                                 chromosome = chr,
                                 start = pos.start,
                                 end = pos.end,
                                 col = "black", fill = "black",
                                 name = "B6 Specific Motif",
                                 shape = "box")
  
  
  # genome and gene track ####
  
  gatrack = GenomeAxisTrack(scale = 0.2, col = "black", fill = "black", lwd = 0.5, fontsize = 6)
  
  
  snp.atrack = AnnotationTrack(nod_snp, genome = "mm10",
                               chromosome = chr,
                               start = pos.start,
                               end = pos.end,
                               stacking = "dense",
                               col = "black",
                               fill = "black",
                               name = "NOD SNPs",
                               collapse = T,
                               shape = "box")
  
  # genes
  grtrack = GeneRegionTrack(range = TxDb.Mmusculus.UCSC.mm10.knownGene,
                            chromosome = chr,
                            # symbol = "Itgav",
                            start = pos.start,
                            end = pos.end,
                            fill = "black",
                            col = "black",
                            fontsize = 4,
                            fontsize.group = 6,
                            fontface = 2,
                            fontcolor = "black",
                            fontcolor.group = "black",
                            name = "UCSC knownGene")
  grtrack
  
  if(length(ranges(grtrack)) > 0 ){
    # convert ensembl transcript id to mgi symbol
    ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", 
                          version = 102, host = "https://nov2020.archive.ensembl.org") 
    
    ranges(grtrack)$symbol = ranges(grtrack)$symbol %>% str_remove(pattern = "\\.[0-9]+")
    id_conversion = getBM(attributes = c('ensembl_transcript_id', 'mgi_symbol'),
                          filters = 'ensembl_transcript_id',
                          values = ranges(grtrack)$symbol, # need to format for use in filters MGI:<ID#>
                          mart = ensembl)
    
    symbols = id_conversion$mgi_symbol
    names(symbols) = id_conversion$ensembl_transcript_id
    
    ranges(grtrack)$symbol = symbols[ranges(grtrack)$symbol]
  }
  
  
  # plot with Highlight ####
  
  # htrack1 <- HighlightTrack(trackList=c(multiome.bw.dtrack
  #                                       ), 
  #                           start = highlight.start,
  #                           end = highlight.end, #83694022-83694322
  #                           chromosome = chr %>% str_remove(pattern = "chr"))
  # htrack2 <- HighlightTrack(trackList = c(gv.bw.dtrack,
  #                                         gv.consensus.atrack), 
  #                           start = highlight.start,
  #                           end = highlight.end, 
  #                           chromosome = chr %>% str_remove(pattern = "chr"))
  
  htrack3 <- HighlightTrack(trackList = c(multiome.bw.dtrack,
                                          egrn.track,
                                          grtrack,
                                          multiome.consensus.atrack,
                                          gv.bw.dtrack,
                                          gv.consensus.atrack,
                                          motif.atrack), 
                            start = highlight.start,
                            end = highlight.end, 
                            col = "grey90",
                            chromosome = chr %>% str_remove(pattern = "chr"))
  
  pdf(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_",regulon,".pdf"),
      width = 4,
      height = 5.5)
  plotTracks(   c(gatrack,
                  htrack3,
                  snp.atrack),
                
                # c(gatrack,
                #            htrack1,
                #            egrn.track,
                #            grtrack,
                #            multiome.consensus.atrack, 
                #            htrack2,
                #            motif.atrack,
                #            snp.atrack),
                
                
                transcriptAnnotation = "symbol", 
                collapseTranscripts = "longest",
                from = pos.start, 
                to = pos.end,
                background.title = "white",
                fontface = 2,
                fontsize = 6,
                #margin = 20,
                #just.group = "below",
                fontcolor = "black",
                col.axis="black")
  dev.off()
  
}