library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)

# mm10 to nod chr conversion

chr_conv = read.delim("/Scottbrowne/seq/tmp/devoes/liftOff/nod-mm10_chromEquiv.csv",
                      header = F, 
                      sep = ",",
                      col.names = c("mm10","NOD"))

# regions in CNV or HetSNP to check where hits are for the ones with thousands
hit_coords = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/Copies_InCNV_or_WithHetSNP/copy_coords.bed",
                    header = F,
                    col.names = c("chr","start","end","peakID","strand","copyNum"))

# idd
idd = read.delim("/Scottbrowne/members/smd/Projects/Idd.bed",
                 header = F,
                 col.names = c("chr","start","end","Idd"))
idd = makeGRangesFromDataFrame(idd, keep.extra.columns = T, starts.in.df.are.0based = T)

# load number of blast hits per peak ####

hits = read.delim("/Scottbrowne/seq/tmp/devoes/SD032/blast_all_peaks/hits_per_peak.csv",
                  header = F,
                  sep = ",")
hits = pivot_wider(hits, names_from = V2 , values_from = V3)
colnames(hits)[1] = "peakID"
hits = hits %>% mutate(diff = NOD - mm10, 
                       chr = str_remove(peakID, "_.*"), 
                       start = str_remove(peakID, ".*_") %>% as.numeric(), 
                       end = start + 301)
hits = makeGRangesFromDataFrame(hits, keep.extra.columns = T, starts.in.df.are.0based = T)

quantile(hits$diff, probs = seq(0,1,0.1))
quantile(hits$NOD, probs = seq(0,1,0.1))
quantile(hits$mm10, probs = seq(0,1,0.1))

hits[hits$peakID == "chr1_171561902"] # cd244
length(which(hits$NOD > hits$mm10))
length(which(hits$NOD < hits$mm10))

ggplot(hits %>% as_tibble(), aes(x = mm10, y = NOD)) + 
  geom_point() + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed',color='red') + 
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) # c(0,100)


# get non-zero differences in blast hits ####
hits_nz = hits[hits$diff !=0]
#hits_nz %>% as_tibble() %>% arrange(diff)

ggplot(hits_nz %>% as_tibble(), aes(x = mm10, y = NOD)) + 
  geom_point() + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed',color='red') + 
  coord_cartesian(xlim = c(0,50), ylim = c(0,50)) # c(0,100)

table(abs(hits_nz$diff) < 100)

# excessive num hits ####

hit_coords[hit_coords$peakID == "chr6_130312322" & grepl("OW", hit_coords$chr), 'chr'] %>% table()
hits[hits$peakID == "chr6_130311867"]
hit_coords[hit_coords$peakID == "chr6_130311867" ,] 


# mhc ####
hits %>% as_tibble() %>% filter(grepl("chr17_352", peakID)) %>% print(n=Inf)

hits_nz %>% as_tibble() %>% filter(abs(diff) > 0 & grepl("chr17",peakID))
hits_nz %>% as_tibble() %>% filter(NOD< 50 & mm10 < 50 & abs(diff) > 0 & grepl("chr17",peakID)) %>% print(n=Inf)

hits_nz[queryHits(findOverlaps(query = hits_nz, 
                               subject = idd[idd$Idd == "Idd1(MHC)"],
                               type = "any",
                               select = "all"))] %>%
  as_tibble() %>% print(n=Inf)

# klra
hits_nz[queryHits(findOverlaps(query = hits_nz, 
                               subject = idd[idd$Idd == "Idd6.AM"],
                               type = "any",
                               select = "all"))] %>%
  as_tibble() %>% print(n=Inf)

# cd244
hits_nz %>% as_tibble() %>% filter(grepl("chr1_1715",peakID))

hits_nz %>% as_tibble() %>% filter(diff > 2 & diff < 15) %>% print(n = Inf)
hits_nz %>% as_tibble() %>% filter(NOD< 50 & mm10 < 50 & abs(diff) > 0) %>% print(n=Inf)

# check only peaks with a HetSNP ####
hSNP = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/ConsensusPeaks_intersect_AnyHetSNP.bed",
                  col.names = c("chr","start","end","peakID"), header = F) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = T)

# MHC
hits_nz[queryHits(findOverlaps(query = hits_nz, 
                               subject = hSNP[queryHits(findOverlaps(query = hSNP, 
                                                                     subject = idd[idd$Idd == "Idd1(MHC)"],
                                                                     type = "any",
                                                                     select = "all"))],
                               type = "any",
                               select = "all"))] %>%
  as_tibble() %>% print(n=Inf)

# klra
hits_nz[queryHits(findOverlaps(query = hits_nz, 
                               subject = hSNP[queryHits(findOverlaps(query = hSNP, 
                                                                     subject = idd[idd$Idd == "Idd6.AM"],
                                                                     type = "any",
                                                                     select = "all"))],
                               type = "any",
                               select = "all"))] %>%
  as_tibble() %>% print(n=Inf)

# any
hits_nz[queryHits(findOverlaps(query = hits_nz, 
                               subject = hSNP,
                               type = "any",
                               select = "all"))] %>%
  as_tibble() %>% arrange(diff) %>% print(n=Inf)

# ####
cnv = read.delim("/Scottbrowne/seq/tmp/devoes/SD029/ConsensusPeaks_intersect_AnyHetSNP_Cahan-Graubert-Alriyami.bed",
                 header = F,
                 col.names = c("chr","start","end","peakID")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T,
                           starts.in.df.are.0based = T)


# mm10.p6 ####
p6.alias = read.delim("/Scottbrowne/seq/tmp/devoes/mm10.p6/mm10.p6.chromAlias.txt")
hits.p6 = read.delim("/Scottbrowne/seq/tmp/devoes/mm10.p6/Copies_InCNV_or_WithHetSNP/copy_coords.tsv",
                     header = F)

hits.p6[hits.p6$V4 == "chr6_130315243",'V1'] %>% table()

p6.alias[grepl("IDD",p6.alias[,2]),]
hits.p6[hits.p6$V4 == "chr6_130312322",'V1'] %>% table() %>% sum()

hits.p6[hits.p6$V4 == "chr6_130312322" & 
          (hits.p6$V1 %in% p6.alias[grepl("IDD",p6.alias[,2]),1]),
        'V1'] %>% table() #  chr6_JH584264_alt IDD6_3

# diff acc peaks intersect ####
dachr = lapply(list.files("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5",
                          pattern = "NOD-B6.*all.bed",
                          full.names = T), 
               rtracklayer::import)
dachr = do.call(c, dachr)
hits_nz[queryHits(findOverlaps(hits_nz, dachr,type="any",select = "all"))] %>% 
  as_tibble() %>%
  unique() %>%
  filter(mm10 > 1 & mm10 < 10 & diff < 0) %>% 
  arrange(peakID) %>%
  print(n=Inf)

# ccl21/ccl27 ####
hits[grepl("chr4_4[1-2][0-9][0-9][0-9][0-9][0-9][0-9]",hits$peakID)] %>% as_tibble() %>% print(n= Inf)
# cd200r1/2/3/4 ####
hits[grepl("chr16_4[4-5][0-9][0-9][0-9][0-9][0-9][0-9]",hits$peakID)] %>% as_tibble() %>% print(n= Inf)
