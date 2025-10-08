library(dplyr)
library(stringr)
library(edgeR)
library(limma)
library(GenomicRanges)


set.seed(75)

# 1. set paths ####

half_width = 150

work_dir = "/Scottbrowne/members/smd/Projects/SD029"

# inputs

metadata_file = paste0(work_dir,"/metadata/low-noise_metadata.tsv")

data_path = "/scratch2/devoes/SD029"
consensus_peak_file = paste0(data_path, # read in RDS to save conversion steps?
                             "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt",
                             "/globallambda/",
                             half_width, 
                             "bp/All_Samples.fwp.filter.non_overlapping.rds")



# outputs

out_data = paste0(work_dir, "/data")

out_figures = paste0(work_dir, "/figures")

# 2. set parameters ####

padj = 0.1

# 3. load data ####

# load sample data

metadata = read.delim(metadata_file)

# set bam paths (NOD and B6 are different)

bam_files = paste0(data_path,
                   "/nfcore_outs/",
                   metadata$Background,
                   "/bowtie2/merged_library/",
                   metadata$nf.core_base,
                   ".sorted.bam")
ind_NOD = metadata$Background == "NOD"
bam_files[ind_NOD] = bam_files[ind_NOD] %>% 
  str_replace(pattern = "bowtie2/merged_library", replacement = "shifted/bowtie2") %>% 
  str_replace(pattern = ".sorted.bam", replacement = ".coordinate-sorted.bam")

# load peak GRanges (1-based)

if(tools::file_ext(consensus_peak_file) == "rds"){
  peak_gr = readRDS(consensus_peak_file)
} else if(tools::file_ext(consensus_peak_file) == "bed") {
  peak_gr = rtracklayer::import(consensus_peak_file)
}
peak_gr$peak_id = paste0(seqnames(peak_gr),"_",start(peak_gr) - 1) # give generic name to peaks

# and convert to saf dataframe
peak_saf = data.frame(peak_gr)
peak_saf = data.frame(
  GeneID = peak_saf$peak_id,
  Chr = peak_saf$seqnames,
  Start = peak_saf$start,
  End = peak_saf$end,
  Strand = peak_saf$strand
)


# 4. count reads in peaks ####

# counts = featureCounts(files = bam_files,
#                        annot.ext = peak_saf, # saf format required
#                        isPairedEnd = T,
#                        countReadPairs = F, # count reads to get biological event of interest (Tn5 insertion)
#                        tmpDir = "/scratch/smd",
#                        allowMultiOverlap = T,
#                        nthreads = 8)
# 

counts = readRDS(paste0(out_data, "/ArchR",half_width,"bpExt_featureCounts.rds"))

# keep only peaks in model

v = readRDS(paste0(out_data,"/v_",half_width,"bp.rds"))
peakIDs = v$genes$GeneID

counts.mtx = counts$counts
counts.mtx = counts.mtx[peakIDs,]

# 5. Calculate size factors for bigwig #####
# https://www.biostars.org/p/413626/#414440

libSize = colSums(counts.mtx)
scaleFactor = 1/(libSize/1e6) # do we want per million reads in peak or per 100k or 10k reads in peaks? to help with visualization? million is pretty standard though
df = data.frame(sample = names(scaleFactor),recipScale = as.vector(scaleFactor))
df = df %>% mutate(sample = str_remove(sample, pattern = ".sorted.bam|.coordinate-sorted.bam"))

write.table(df, file = "/Scottbrowne/members/smd/Projects/SD029/sample_lists/bw-customScale_input.csv",
            quote = F, col.names = F, row.names = F,
            sep = ",")
