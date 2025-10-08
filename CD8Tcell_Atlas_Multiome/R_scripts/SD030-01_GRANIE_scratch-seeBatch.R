library(GRaNIE)
library(stringr)
library(dplyr)
library(ggplot2)

# paths and parameters ####
half_width = 150
atac_cpm_filt = 3
cpm_filt_rna = 1
num_groups_cpm_filt = 1
genomeAssembly = "mm10"

data_path = "/Scottbrowne/seq/tmp"
proj_path = "/Scottbrowne/members/smd/Projects/SD030"

meta_data_path = paste0(proj_path, "/metadata/low-noise_metadata.txt")
meta_data = read.delim(meta_data_path, 
                       header = T,
                       sep = "\t")
table(meta_data$archr_peak_group)

rna_mtx_path = paste0(data_path, "/SD030/nfcore_outs/rna/",
                      "star_salmon/salmon.merged.gene_counts.tsv") # salmon.merge.transcript_counts.tsv
tx2gene = read.delim(paste0(data_path,"/SD030/nfcore_outs/rna/star_salmon/tx2gene.tsv"),
                     header = F)
bam_path = meta_data$nfcore_bam_path
peak_path = paste0(data_path, "/devoes/SD030/macs2_bed",
                  "/archr-iterative-merge_peaks-by-group_variableExt/globallambda/",
                  half_width,
                  "bp/All_Samples.fwp.filter.non_overlapping.renamed.bed"
                  )

# filter RNA ####

rna_counts = read.delim(rna_mtx_path)
rna_counts = rna_counts[,grepl(pattern = paste0(collapse = "|", c('gene', paste0(meta_data$SampleID, "_RNA"))), 
                               colnames(rna_counts))]
rownames(rna_counts) = rna_counts$gene_id
colnames(rna_counts) = colnames(rna_counts) %>% str_remove(pattern = "_RNA")
# determine cpm approx 10 counts
rna_lib_size = colSums(rna_counts[,!grepl("gene",colnames(rna_counts))])
rna_cpm = rna_counts[,!grepl("gene",colnames(rna_counts))] / rna_lib_size * 1e6
colnames(rna_cpm) = colnames(rna_cpm) %>% str_remove(pattern = "_RNA")
quantile(rna_lib_size, probs = seq(0,1,0.1)) # what do we want to do about the widely varied lib sizes? Do we want to use a flat counts cutoff instead?
# get avg cpm by group
if(all(colnames(rna_cpm) == meta_data$sampleID)){
  rna_cpm_group = data.frame(t(rna_cpm), group = meta_data$archr_peak_group)
}
all(meta_data$SampleID == rownames(rna_cpm_group))
rna_cpm_group = rna_cpm_group %>% group_by(group) %>% summarise(across(.cols = everything(), mean))
# check distribution and get passing peaks
rm_group_rna = !grepl(pattern = "group",colnames(rna_cpm_group))
hist(log10(rna_cpm_group[,rm_group_rna] %>%
       as.vector() %>% unlist() + 0.5),
     breaks = 200)
quantile(rna_cpm_group[,rm_group_rna] %>%
           as.vector() %>% unlist(),
         probs=seq(0,1,0.05))
# want gene to be prominent in at least one group
keep_genes = colSums(rna_cpm_group[,rm_group_rna] > cpm_filt_rna) >= num_groups_cpm_filt

# or do we want to do flat 10 counts?
#rowSums(rna_counts[,!grepl('gene',colnames(rna_counts))] > 10 ) > 3

table(keep_genes)
keep_genes = names(keep_genes)[keep_genes]

rna_counts_filt = rna_counts[keep_genes,!grepl('gene',colnames(rna_counts))]

# format for GRANIE
rna_counts_filt = as_tibble(rna_counts_filt) %>%
  mutate(ENSEMBL = keep_genes) %>%
  relocate(ENSEMBL)
colnames(rna_counts_filt)

?limma::removeBatchEffect()

# filter ATAC ####

peak_gr = rtracklayer::import(peak_path)

# check if non-standard chr present
seqnames(peak_gr) %>% table()

# peaks to saf
peak_saf = data.frame(peak_gr)
peak_saf = data.frame(
  GeneID = peak_saf$name,
  Chr = peak_saf$seqnames,
  Start = peak_saf$start, # start is already 1-based from granges import
  End = peak_saf$end,
  Strand = peak_saf$strand
)


# get counts per peak

# atac_counts = Rsubread::featureCounts(files = bam_path,
#                        annot.ext = peak_saf, # saf format required
#                        isPairedEnd = T,
#                        countReadPairs = F, # count reads to get biological event of interest (Tn5 insertion)
#                        tmpDir = "/scratch2/devoes/tmp",
#                        allowMultiOverlap = T,
#                        nthreads = 8)
# saveRDS(atac_counts, file = paste0(proj_path, "/data/ArchR",half_width, "bpExt_featureCounts.rds"))


atac_counts = readRDS(paste0(proj_path, "/data/ArchR",half_width,"bpExt_featureCounts.rds"))
atac_counts_mtx = atac_counts$counts
colnames(atac_counts_mtx) = colnames(atac_counts_mtx) %>% str_remove(pattern = "_ATAC_.*")

# get cpm
atac_lib_size = colSums(atac_counts_mtx)
atac_cpm = atac_counts_mtx / atac_lib_size * 1e6
quantile(atac_lib_size, probs = seq(0,1,0.1)) # what do we want to do about the widely varied lib sizes? Do we want to use a flat counts cutoff instead?
# get avg cpm by group
atac_cpm_group = data.frame(t(atac_cpm), group = meta_data$archr_peak_group)
all(meta_data$SampleID == rownames(atac_cpm_group))
atac_cpm_group = atac_cpm_group %>% group_by(group) %>% summarise(across(.cols = everything(), mean))
# check distribution and get passing peaks
rm_group = !grepl(pattern = "group",colnames(atac_cpm_group))
hist(atac_cpm_group[,rm_group] %>%
       as.vector() %>% unlist(),
     breaks = 200)
quantile(atac_cpm_group[,rm_group] %>%
           as.vector() %>% unlist(),
         probs=seq(0,1,0.05))
# want peak to be prominent in at least one group
keep_peaks = colSums(atac_cpm_group[,rm_group] > atac_cpm_filt) >= num_groups_cpm_filt
table(keep_peaks)
keep_peaks = names(keep_peaks)[keep_peaks]

# or do we want to do flat 10 counts?
#table(rowSums(atac_counts_mtx > 10 ) > 3)


# filter by cpm
atac_counts_mtx_filt = atac_counts_mtx[keep_peaks,]

# format for GRANIE
atac_counts_filt = as_tibble(atac_counts_mtx_filt) %>%
  mutate(peakID = keep_peaks) %>%
  relocate(peakID)
colnames(atac_counts_filt)

peakIDs = paste0(peak_saf$Chr,
                 ":",
                 peak_saf$Start-1, #adjust to zero based (bed)
                 "-",
                 peak_saf$End)
names(peakIDs) = peak_saf$GeneID
atac_counts_filt$peakID = peakIDs[atac_counts_filt$peakID]


# initialize GRANIE objet #####

objectMetadata.l = list(name = paste0("CD8 T cells in Cancer and Chronic Infection"),
                        file_peaks = peak_path,
                        file_rna = rna_mtx_path, 
                        file_sampleMetadata = meta_data_path,
                        genomeAssembly = genomeAssembly)
GRN = initializeGRN(objectMetadata = objectMetadata.l,
                    outputFolder = paste0(proj_path, "/data/GRANIE"),
                    genomeAssembly = genomeAssembly)

GRN

GRN = addData(GRN, 
              counts_peaks = atac_counts_filt,
              normalization_peaks = "limma_quantile", # using limma quantile to match our GV workflow
              idColumn_peaks = "peakID",
              counts_rna = rna_counts_filt, 
              normalization_rna = "limma_quantile",
              idColumn_RNA = "ENSEMBL",
              genomeAnnotationSource = "AnnotationHub",
              sampleMetadata = meta_data, # check that _ATAC and _RNA are removed from count mtxs
              forceRerun = TRUE)

GRN

GRN = plotPCA_all(GRN, data = c("rna","peaks"), topn = 500, type = "normalized", plotAsPDF = T,
                  #pages = c(2, 3, 14),
                  forceRerun = TRUE)

# distance/similarity metrics ####

# RNA MDS
rna_norm_counts = getCounts(GRN, type = "rna", includeIDColumn = F)
rna_eucl.dist = dist(t(rna_norm_counts)) # need logtransformed, normalized counts
rna_mds_fit = cmdscale(rna_eucl.dist, eig = T, k =2)
rna_mds = data.frame(rna_mds_fit$points)
rna_mds$SampleID = rownames(rna_mds) 
rna_mds = left_join(rna_mds, meta_data, by = "SampleID")
ggplot(rna_mds, aes(x = X1, y = X2, shape = Phenotype, color = Publication)) +
  geom_point(size = 4) +
  theme_classic()

# ATAC MDS
atac_norm_counts = getCounts(GRN, type = "peaks", includeIDColumn = F) 
atac_eucl.dist = dist(t(atac_norm_counts)) # need logtransformed, normalized counts
atac_mds_fit = cmdscale(atac_eucl.dist, eig = T, k =2)
atac_mds = data.frame(atac_mds_fit$points)
atac_mds$SampleID = rownames(atac_mds) 
atac_mds = left_join(atac_mds, meta_data, by = "SampleID")
ggplot(atac_mds, aes(x = X1, y = X2, shape = Phenotype, color = Publication)) +
  geom_point(size = 4) +
  theme_classic()
