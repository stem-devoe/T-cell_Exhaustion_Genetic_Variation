library(dplyr)
library(stringr)
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)


# parameters ####

regions = read.delim("/Scottbrowne/members/smd/Projects/SD029/plotRegions/QC_regions.txt")
regions$chr = str_remove(regions$IGV, pattern = ":.*")
regions$start = str_remove(regions$IGV, pattern = ".*:") %>% 
  str_remove(pattern = "-.*") %>% 
  str_remove_all(pattern = ",") %>% 
  as.numeric()
regions$end = str_remove(regions$IGV, pattern = ".*-") %>% 
  str_remove_all(pattern = ",") %>% 
  as.numeric()

for(i in seq_len(nrow(regions))){
  genome = "mm10"
  chr = regions$chr[i] # "chr2"
  gene = regions$gene_locus[i] # "Itgav"
  pos.start = regions$start[i] #83693000  
  pos.end = regions$end[i] # 83750000 
  gv.bw.lim = c(0,200)
  highlight.start = regions$highlight_start[i]	# 83694022
  highlight.end = regions$highlight_end[i]
  
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
  gv_bw_files = list.files("/scratch2/devoes/nfcore_GV/bigwig_frip_groupAvg",
                           full.names = T,
                           pattern = ".bw") %>% grep(pattern = "Naive|CART", value = T)
  gv_bw_files = c(gv_bw_files[grepl("CART", gv_bw_files)], gv_bw_files[grepl("Naive", gv_bw_files)])
  names(gv_bw_files) = gv_bw_files %>% basename() %>% str_remove(pattern = "_frip.*")
  
  
  # peaks
  gv_consensus = rtracklayer::import("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.renamed.bed")
  gv_consensus = gv_consensus[subjectHits(findOverlaps(GRanges(seqnames = chr,
                                                               ranges = IRanges(start = pos.start, end = pos.end),
                                                               strand = "*"),
                                                       gv_consensus)),]
  
  
  # NOD v B6 tracks ####
  
  # bigwigs
  
  gvDataTrack = function(bw_file, fill_colors, chr, bw.lim){
    DataTrack(bw_file, genome = "mm10",
              type = "h",
              window  = -1, # add smoothing effect
              chromosome = chr,
              ylim = bw.lim,
              yTicksAt = bw.lim[2],
              fontsize.legend = 6,
              col = fill.colors[bw_file %>% basename() %>% str_remove(pattern = "_frip.*")])
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
                                        shape = "box")
  
  
  # genome and gene track ####
  
  gatrack = GenomeAxisTrack(scale = 0.2, col = "black", fill = "black", lwd = 0.5, fontsize = 6)
  
  
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
  #grtrack
  
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
  
  # plot ####
  
  pdf(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_",gene,".pdf"),
      width = 3,
      height = 2)
  plotTracks(c(gatrack,
               gv.consensus.atrack,
               gv.bw.dtrack,
               grtrack), #snp.atrack
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
  
  
  # plot with Highlight ####
  
  if(gene == "Pdcd1"){
    htrack1 <- HighlightTrack(trackList=c(gv.bw.dtrack), 
                              start = highlight.start,
                              end = highlight.end, 
                              chromosome = chr %>% str_remove(pattern = "chr"))
    pdf(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_",gene,"_highlight.pdf"),
        width = 3,
        height = 2)
    plotTracks(c(gatrack,
                 gv.consensus.atrack,
                 htrack1,
                 grtrack),
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

  
}
