library(dplyr)
library(GenomicRanges)

# motif distance

lost_peaks = rtracklayer::import("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5/NOD-B6_CART.down.bed")

# load motif of interest
moi = "SP1"
moi.files = list.files(path = "/scratch2/devoes/SD029/HOCOMOCO/categorized/",
                       pattern = paste0("ConsensusPeaks_HOCOMOCO_", moi,".*.categorized.bed"),
                       full.names = T)
moi.grl = lapply(moi.files, read.delim, header = F, col.names = c("chr","start","stop","ID","strand","posseq","specificity")) 
moi.grl = lapply(moi.grl, makeGRangesFromDataFrame, keep.extra.columns = T,starts.in.df.are.0based = T)
moi.gr = purrr::reduce(moi.grl, c)
moi.gr = moi.gr[!grepl("NOD", moi.gr$specificity)]

moi.gr.b6 = moi.gr[grepl("B6", moi.gr$specificity)]
moi.gr.b6.lost = moi.gr[queryHits(findOverlaps(query = moi.gr.b6, subject = lost_peaks))]


neighbor = "NR4"
neighbor.files = list.files(path = "/scratch2/devoes/SD029/HOCOMOCO/categorized/",
                            pattern = paste0("ConsensusPeaks_HOCOMOCO_", neighbor,".*.categorized.bed"),
                            full.names = T)
neighbor.grl = lapply(neighbor.files, read.delim, header = F, col.names = c("chr","start","stop","ID","strand","posseq","specificity")) 
neighbor.grl = lapply(neighbor.grl, makeGRangesFromDataFrame, keep.extra.columns = T,starts.in.df.are.0based = T)
neighbor.gr = purrr::reduce(neighbor.grl,c)
neighbor.gr = neighbor.gr[!grepl("NOD", neighbor.gr$specificity)]

# get motif midpoint

nearest.all = distanceToNearest(moi.gr, neighbor.gr)
quantile(nearest.all@elementMetadata$distance)
max(nearest.all@elementMetadata$distance)
#plot(density(nearest.all@elementMetadata$distance))
table(nearest.all@elementMetadata$distance < 300)

nearest.b6 = distanceToNearest(moi.gr.b6, neighbor.gr)
quantile(nearest.b6@elementMetadata$distance)
max(nearest.b6@elementMetadata$distance)
#plot(density(nearest.b6@elementMetadata$distance))
table(nearest.b6@elementMetadata$distance < 300)

nearest.b6.lost = distanceToNearest(moi.gr.b6.lost, neighbor.gr)
quantile(nearest.b6.lost@elementMetadata$distance)
max(nearest.b6.lost@elementMetadata$distance)
table(nearest.b6.lost@elementMetadata$distance < 300)
  
# ~20% within peak width of RUNX
# ~50% within peak width of ETS