library(dplyr)
library(motifmatchr)
library(memes)
library(GenomicRanges)
library(stringr)
library(TFBSTools)

# adjust options so results are NOT written in scientific notation ####
if(options()$scipen != 9999){
  options(scipen = 9999)
}

# motif formatting ####

motif_file = "/Scottbrowne/members/smd/motifs/HOCOMOCO/v11/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"
# read meme to universal motif format - we want to run one motif at a time
motif.univ = universalmotif::read_meme(motif_file)
motif.names = lapply(motif.univ, FUN = function(x){x@name}) %>% unlist()
names(motif.univ) = motif.names

# 1. Paths ####

proj_dir = "/Scottbrowne/members/smd/Projects/SD029/"
scratch_dir = "/scratch2/devoes/SD029"

if(!dir.exists(paste0(scratch_dir,"/HOCOMOCO/fimp_outs"))){
  dir.create(paste0(scratch_dir,"/HOCOMOCO/fimo_outs"), recursive =  T)
  dir.create(paste0(scratch_dir,"/HOCOMOCO/runFimo"), recursive =  T)
  dir.create(paste0(scratch_dir,"/HOCOMOCO/categorized"), recursive =  T)
}

fa_files = list.files(paste0(scratch_dir,
                             "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp"),
                      full.names = T,
                      pattern = ".fa"
)
peak_files = list.files(paste0(scratch_dir,
                               "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp"),
                        full.names = T,
                        pattern = ".bed") %>% grep(pattern = "renamed", value = T)
# check if peak width is maintained in NOD genome
nod_gr = rtracklayer::import(peak_files %>% grep(pattern = "NOD", value = T))
quantile(width(nod_gr))


meme_path = "/home/devoes/miniconda3/envs/memeSuite/bin"


# 3. find motifs with memesuite FIMO ####

saveMotifSites <- function(motif_name, path2fasta, motifSet, meme_path, fimo_outdir, outdir, outbase, thresh = 1e-04, quiet = F){ 
  
  if(!quiet){
    print(paste0("Running Fimo ", motif_name))
  }
  fimo_res = findMotifSites(motif_name, path2fasta, motifSet, meme_path, fimo_outdir, thresh = thresh)
  if(!quiet){
    print("Saving results")
  }
  saveFimoRDS(motif_name, fimo_res, outdir, outbase)
  if(!quiet){
    print("Converting results to bed format")
  }
  fimo_bed = format2bed(fimo_res)
  if(!quiet){
    print("Saving bed")
  }
  saveMotifBed(fimo_bed, motif_name, outdir,outbase)
  return()
}

findMotifSites <- function(motif_name, path2fasta, motifSet, meme_path, fimo_outdir, thresh){
  fimo_res <- runFimo(
    sequences = path2fasta, # fasta
    motifs = motifSet[motif_name], # can do more than one but run 1 at a time for saving files
    thresh = thresh, # 1e-04 is default
    bfile = "motif",
    outdir = fimo_outdir, #
    parse_genomic_coord = T,
    skip_matched_sequence = F,
    max_strand = F, # want all matches within peak if more than 1
    text = T,
    meme_path = meme_path,
    silent = T)
  return(fimo_res)
}

saveFimoRDS <- function(motif_name, fimo_res, outdir, outbase){
  # remove colons and commna from motif name
  motif_name = motif_name %>% 
    str_replace_all(pattern = "\\?|\\+|:|,|\\||/|\\\\", replacement = ".") 
  outfile = paste0(outdir,"/",outbase,"_",motif_name,".rds")
  saveRDS(fimo_res,outfile)
  return()
}

format2bed <- function(fimo_res){
  
  if(is.null(fimo_res)){
    fimo_bed = data.frame()
   #colnames(fimo_bed) = c('seqnames','start','end','ID')
  }else{
    fimo_bed = data.frame(fimo_res)
    fimo_bed$start = fimo_bed$start - 1 # adjust back to 0-base
    fimo_bed$ID = with(fimo_bed, 
                       paste0(seqnames, ":", start, "_", motif_id, "_",
                              strand, "_", score, "_", matched_sequence))
    fimo_bed = fimo_bed[,c('seqnames','start','end','ID')]
  }

  return(fimo_bed)
}

saveMotifBed <- function(fimo_bed, motif_name, outdir, outbase){
  # remove colons and commna from motif name
  motif_name = motif_name %>% 
    str_replace_all(pattern = "\\?|\\+|:|,|\\||/|\\\\", replacement = ".") 
  outfile = paste0(outdir, "/",outbase,"_",motif_name,".bed")
  write.table(fimo_bed,
              sep = "\t",
              quote = F, col.names = F, row.names = F,
              file = outfile)
  return()
}

## a. B6 sequence ####

lapply(motif.names, FUN = saveMotifSites,
       path2fasta = grep("NOD", fa_files, value = T, invert = T),
       motifSet = motif.univ,
       meme_path = meme_path,
       fimo_outdir = paste0(scratch_dir, "/HOCOMOCO/runfimo"),
       outdir = paste0(scratch_dir, "/HOCOMOCO/fimo_outs"),
       outbase = "ConsensusPeaks_B6_HOCOMOCO")

## b. NOD sequence ####

lapply(motif.names, FUN = saveMotifSites,
       path2fasta = grep("NOD", fa_files, value = T, invert = F),
       motifSet = motif.univ,
       meme_path = meme_path,
       fimo_outdir = paste0(scratch_dir, "/HOCOMOCO/runFimo"),
       outdir = paste0(scratch_dir, "/HOCOMOCO/fimo_outs"),
       outbase = "ConsensusPeaks_NOD_HOCOMOCO",
       quiet = T)



