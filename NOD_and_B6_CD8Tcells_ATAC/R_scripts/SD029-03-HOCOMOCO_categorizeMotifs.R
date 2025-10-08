library(dplyr)
library(stringr)
library(GenomicRanges)

# 0.  ####
scratch_dir = "/scratch2/devoes/SD029"

b6.motif.files = list.files(path = paste0(scratch_dir, "/HOCOMOCO/fimo_outs"),
                            pattern = "ConsensusPeaks_B6_HOCOMOCO_.*.bed",
                            full.names = T)
names(b6.motif.files) = str_match(b6.motif.files, pattern = "HOCOMOCO_(.*).bed") %>% .[,2]

nod.motif.files = list.files(path = paste0(scratch_dir,"/HOCOMOCO/fimo_outs"),
                             pattern = "shifted_from_NOD_SHILTJ.bed",
                             full.names = T)
names(nod.motif.files) = str_match(nod.motif.files, pattern = "HOCOMOCO_(.*)_shifted") %>% .[,2]

moi = intersect(names(nod.motif.files), names(b6.motif.files)) # some files will be empty and I need to account for that


lapply(moi,
       FUN = categorizeMotifs,
       b6.files = b6.motif.files, 
       nod.files = nod.motif.files, 
       outdir = "/scratch2/devoes/SD029/HOCOMOCO/categorized",
       outbase = "ConsensusPeaks_HOCOMOCO")


# FUNCTIONS #####
categorizeMotifs <- function(motif_name, b6.files, nod.files, outdir, outbase){
  
  
  # 0. check for no hits being found ####
  b6.empty = file.size(b6.files[motif_name]) == 0L
  nod.empty = file.size(nod.files[motif_name]) == 0L
  
  if(!b6.empty & !nod.empty){
    print(paste0("Hits in both strains: ", motif_name))
    hitsInBoth(motif_name, b6.files, nod.files, outdir, outbase)
  } else if(!b6.empty){
    print(paste0("Hits only in B6: ", motif_name))
    hitsInOne(motif_name, strain = "B6",
              strain.files = b6.files,
              outdir, outbase)
  } else if (!nod.empty){
    print(paste0("Hits only in NOD: ", motif_name))
    hitsInOne(motif_name, strain = "NOD", 
              strain.files = nod.files,
              outdir, outbase)
  } else{ # no hits in either strain
    print(paste0("No hits in either strain: ", motif_name))
    return()
  }
  
  return()
}

findMotifSignature <- function(strain, gr){
  # signature is: chr:start:strand(?):sequence - should be able to handle strands with _(\\+|-)_ vs motif with a -
  
  if(strain == "B6"){
    b6.info = str_match(gr$name, pattern = "(chr.*:[0-9]*)_.*_(\\+|-)_.*_([AGCTacgt]*)") %>% data.frame()
    b6.signat = with(b6.info, paste(X2, X3, X4,sep = ";")) # use ; to make it easier to parse
    return(b6.signat)
  } else{
    nod.info = str_match(gr$name, pattern = "_(\\+|-)_.*_([AGCTacgt]*)") %>% data.frame()
    # need to get B6 coordinates from gr not from name field
    nod.signat = paste0(seqnames(gr),":",start(gr)-1, ";", nod.info$X2, ";", nod.info$X3) 
    return(nod.signat)
  }
  
  
}

saveBED <- function(categorized_gr, motif_name, outdir, outbase){
  bed = data.frame(categorized_gr)
  bed$start = bed$start - 1 # adjust from 1-base to 0-base
  bed$strand = str_replace(bed$strand, pattern = "\\*", replacement = ".")
  bed = bed[, c('seqnames','start','end','name','strand','signature','category')]
  
  write.table(bed,
              file =  paste0(outdir, "/", outbase, "_", motif_name, "_categorized.bed"),
              sep = "\t",
              quote = F,
              col.names = F,
              row.names = F)
}

hitsInBoth <- function(motif_name, b6.files, nod.files, outdir, outbase){
  # 1. load as gr object ####
  b6_gr = rtracklayer::import(b6.files[motif_name])
  nod_gr = rtracklayer::import(nod.motif.files[motif_name])
  
  # 2. get motif signature ####
  
  b6_gr$signature = findMotifSignature("B6", b6_gr)
  nod_gr$signature = findMotifSignature("NOD", nod_gr)
  
  # 3. find identical motifs between B6 and NOD ####
  
  # and set aside in a unique list
  b6.same = b6_gr$signature %in% nod_gr$signature
  nod.same = nod_gr$signature %in% b6_gr$signature
  
  common_gr = b6_gr[b6.same]
  common_gr$category = "Common"
  
  # 4. find common but polymorphic #####
  
  b6_diff = b6_gr[!b6.same]
  #print(b6_diff)
  nod_diff = nod_gr[!nod.same]
  #print(nod_diff)
  matches = findOverlaps(b6_diff, nod_diff, type = "any", select = "all")
  #unique(queryHits(matches)) %>% length() %>% print()
  #unique(subjectHits(matches)) %>% length() %>% print()
  
  if(length(queryHits(matches)) == 0){
    polymorphic_gr = NULL
  } else{
    polymorphic_gr = b6_diff[queryHits(matches)]
    #print(polymorphic_gr)
    polymorphic_gr$name = paste0(b6_diff$name[queryHits(matches)],";",nod_diff$name[subjectHits(matches)])
    polymorphic_gr$signatureure = paste0(b6_diff$signatureure[queryHits(matches)],";",nod_diff$signatureure[subjectHits(matches)])
    polymorphic_gr$category = "Common_Polymorphic"
  }
  
  
  # 5. strain specific ####
  
  if(length(queryHits(matches)) == 0){
    b6_specific_gr = b6_diff
    b6_specific_gr$category = "B6_Specific"
    nod_specific_gr = nod_diff
    nod_specific_gr$category = "NOD_Specific"
  } else{
    b6_specific_gr = b6_diff[-queryHits(matches)]
    b6_specific_gr$category = "B6_Specific"
    nod_specific_gr = nod_diff[-subjectHits(matches)]
    nod_specific_gr$category = "NOD_Specific"
  }
  
  
  # 6. Regroup and save ####
  
  categorized_gr = c(common_gr, polymorphic_gr, b6_specific_gr, nod_specific_gr)
  
  saveRDS(categorized_gr, paste0(outdir, "/",outbase, "_",motif_name,"_categorized.rds"))
  
  saveBED(categorized_gr, motif_name, outdir, outbase)
  
  return()
}

hitsInOne <- function(motif_name, strain, strain.files, outdir, outbase){
  # 1. load as gr object ####
  
  strain_gr = rtracklayer::import(strain.files[motif_name])
  
  # 2. get motif signature ####
  
  strain_gr$signature = findMotifSignature(strain, strain_gr)
  
  # 3. strain specific ####
  
  strain_gr$category = paste0(strain, "_Specific")
  
  # 6. Regroup and save ####
  
  saveRDS(strain_gr, paste0(outdir, "/",outbase, "_",motif_name,"_categorized.rds"))
  
  saveBED(strain_gr, motif_name, outdir, outbase)
  
  return()
}

