library(stringr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(limma) 
library(GRaNIE)

#NOTE SB2019 PD1loTIM3lo are TUMOR-NONREACTIVE not naive(mislabeled in sample ID but corrected for group and phenotype)

# paths and parameters ####
half_width = 150
atac_cpm_filt = 3
cpm_filt_rna = 0.5
num_groups_cpm_filt = 1
genomeAssembly = "mm10"

data_path = "/Scottbrowne/seq/tmp"
proj_path = "/Scottbrowne/members/smd/Projects/SD030"

sample_data_path = paste0(proj_path, "/metadata/granie_metadata.txt") # omits ScottBrowne2016_Terminal_Infection_None group due to 1 sample
sample_data = read.delim(sample_data_path, 
                         header = T,
                         sep = "\t")
table(sample_data$archr_peak_group)

rna_mtx_path = paste0(data_path, "/SD030/nfcore_outs/rna/",
                      "star_salmon/salmon.merged.gene_counts.tsv") # salmon.merge.transcript_counts.tsv
tx2gene = read.delim(paste0(data_path,"/SD030/nfcore_outs/rna/star_salmon/tx2gene.tsv"),
                     header = F)
bam_path = sample_data$nfcore_bam_path
peak_path = paste0(data_path, "/devoes/SD030/macs2_bed",
                   "/archr-iterative-merge_peaks-by-group_variableExt/globallambda/",
                   half_width,
                   "bp/All_Samples.fwp.filter.non_overlapping.renamed.bed"
)

# prepare RNA ####

# load counts
rna_counts = read.delim(rna_mtx_path)
colnames(rna_counts) = colnames(rna_counts) %>% str_remove(pattern = "_RNA")
rna_counts = rna_counts[,c('gene_id',sample_data$SampleID)] # remove extra samples 
rownames(rna_counts) = rna_counts$gene_id

# determine cpm for approx 10 counts
rna_lib_size = colSums(rna_counts[,!grepl("gene",colnames(rna_counts))])
quantile(rna_lib_size, probs = seq(0,1,0.1)) # what do we want to do about the widely varied lib sizes? Do we want to use a flat counts cutoff instead?

rna_cpm = edgeR::cpm(rna_counts[,!grepl("gene",colnames(rna_counts))], log = FALSE) 

# get avg cpm by group
if(all(colnames(rna_cpm) == sample_data$sampleID)){
  rna_cpm_group = data.frame(t(rna_cpm), group = sample_data$archr_peak_group)
} else{
  print("fix")
}
all(sample_data$SampleID == rownames(rna_cpm_group))
rna_cpm_group = rna_cpm_group %>% group_by(group) %>% summarise(across(.cols = everything(), mean))
# check distribution and get passing peaks
rm_group_rna = !grepl(pattern = "group",colnames(rna_cpm_group))
hist(log10(rna_cpm_group[,rm_group_rna] %>%
             as.vector() %>% unlist() + 0.01),
     breaks = 200)
quantile(rna_cpm_group[,rm_group_rna] %>%
           as.vector() %>% unlist(),
         probs=seq(0,1,0.05))
# want gene to be prominent in at least one group
keep_genes = colSums(rna_cpm_group[,rm_group_rna] > cpm_filt_rna) >= num_groups_cpm_filt
table(keep_genes) # 19281

# or do we want to do flat 10 counts?
#table(rowSums(rna_counts[,!grepl('gene',colnames(rna_counts))] > 10 ) > 3) # 17293

keep_genes = names(keep_genes)[keep_genes]

rna_counts_filt = rna_counts[keep_genes,!grepl('gene',colnames(rna_counts))]

# normalize
model_rna = edgeR::DGEList(counts = rna_counts_filt,
                           group = sample_data$archr_peak_group,
                           genes = rownames(rna_counts_filt),
                           samples = sample_data) 
design_rna <- model.matrix(~archr_peak_group + 0, # set no intercept
                           # covariate , # do we want to add covariates? tumor?
                           model_rna$samples) 
voom_rna = voom(model_rna, design_rna, plot = T, normalize = "quantile")

# correct for batch effects
rna_norm_counts = voom_rna$E
rna_granie_mtx = removeBatchEffect(voom_rna$E, batch = model_rna$samples$Publication)

rna_eucl.dist = dist(t(rna_granie_mtx)) # need logtransformed, normalized counts
rna_mds_fit = cmdscale(rna_eucl.dist, eig = T, k =2)
rna_mds = data.frame(rna_mds_fit$points)
rna_mds$SampleID = rownames(rna_mds) 
rna_mds = left_join(rna_mds, sample_data, by = "SampleID")
ggplot(rna_mds, aes(x = X1, y = X2, shape = Phenotype, color = Publication)) +
  geom_point(size = 4) +
  theme_classic()

# format for GRANIE
rna_granie_mtx = as_tibble(rna_granie_mtx) %>%
  mutate(ENSEMBL = rownames(rna_granie_mtx)) %>%
  relocate(ENSEMBL)
colnames(rna_granie_mtx)


# prepare ATAC ####

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
#                                       annot.ext = peak_saf, # saf format required
#                                       isPairedEnd = T,
#                                       countReadPairs = F, # count reads to get biological event of interest (Tn5 insertion)
#                                       tmpDir = "/Scottbrowne/seq/tmp/devoes/tmp",
#                                       allowMultiOverlap = T,
#                                       nthreads = 8)
# saveRDS(atac_counts, file = paste0(proj_path, "/data/ArchR",half_width, "bpExt_featureCounts.rds"))


atac_counts = readRDS(paste0(proj_path, "/data/ArchR",half_width,"bpExt_featureCounts.rds"))
atac_counts_mtx = atac_counts$counts
colnames(atac_counts_mtx) = colnames(atac_counts_mtx) %>% str_remove(pattern = "_ATAC_.*")

# get cpm
atac_lib_size = colSums(atac_counts_mtx)
atac_cpm = atac_counts_mtx / atac_lib_size * 1e6
quantile(atac_lib_size, probs = seq(0,1,0.1)) # what do we want to do about the widely varied lib sizes? Do we want to use a flat counts cutoff instead?
# get avg cpm by group
atac_cpm_group = data.frame(t(atac_cpm), group = sample_data$archr_peak_group)
all(sample_data$SampleID == rownames(atac_cpm_group))
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
atac_counts_filt = atac_counts_mtx[keep_peaks,]

# normalize
model_atac = edgeR::DGEList(counts = atac_counts_filt,
                            group = sample_data$archr_peak_group,
                            genes = rownames(atac_counts_filt),
                            samples = sample_data) 
design_atac <- model.matrix(~archr_peak_group + 0, # set no intercept
                            # covariate , # do we want to add covariates? tumor?
                            model_atac$samples) 
voom_atac = voom(model_atac, design_atac, plot = T, normalize = "quantile")

# correct for batch effects
atac_norm_counts = voom_atac$E
atac_granie_mtx = removeBatchEffect(voom_atac$E, batch = model_atac$samples$Publication)

atac_eucl.dist = dist(t(atac_granie_mtx)) # need logtransformed, normalized counts
atac_mds_fit = cmdscale(atac_eucl.dist, eig = T, k =2)
atac_mds = data.frame(atac_mds_fit$points)
atac_mds$SampleID = rownames(atac_mds) 
atac_mds = left_join(atac_mds, sample_data, by = "SampleID")
ggplot(atac_mds, aes(x = X1, y = X2, shape = Phenotype, color = Publication)) +
  geom_point(size = 4) +
  theme_classic()


# format for GRANIE
atac_granie_mtx = as_tibble(atac_granie_mtx) %>%
  mutate(peakID = rownames(atac_granie_mtx)) %>%
  relocate(peakID)
colnames(atac_granie_mtx)

peakIDs = paste0(peak_saf$Chr, # need chr:start-end
                 ":",
                 peak_saf$Start-1, #adjust to zero based (bed)
                 "-",
                 peak_saf$End)
names(peakIDs) = peak_saf$GeneID
atac_granie_mtx$peakID = peakIDs[atac_granie_mtx$peakID]

# initialize GRANIE objet #####

objectMetadata.l = list(name = paste0("CD8 T cells in Cancer and Chronic Infection"),
                        file_peaks = peak_path,
                        file_rna = rna_mtx_path, 
                        file_sampleMetadata = sample_data_path,
                        genomeAssembly = genomeAssembly)
GRN = initializeGRN(objectMetadata = objectMetadata.l,
                            outputFolder = paste0(proj_path, "/data/GRANIE"),
                            genomeAssembly = genomeAssembly)

GRN

GRN = addData(GRN, 
                      counts_peaks = atac_granie_mtx,
                      normalization_peaks = "none", # normalization done before to handle batch effects
                      idColumn_peaks = "peakID",
                      counts_rna = rna_granie_mtx, 
                      normalization_rna = "none", # normalization done before to handle batch effects
                      idColumn_RNA = "ENSEMBL",
                      #EnsemblVersion = 102, ##
                      #genomeAnnotationSource = "biomaRt",
                      genomeAnnotationSource = "AnnotationHub", # set annotation to mm39 even though I have mm10 in 
                      sampleMetadata = sample_data, # check that _ATAC and _RNA are removed from count mtxs
                      forceRerun = TRUE)


GRN
GRN@annotation$genes[GRN@annotation$genes$gene.name == "Itgav" & !is.na(GRN@annotation$genes$gene.name),]
GRN@annotation$genes[GRN@annotation$genes$gene.name == "Icos" & !is.na(GRN@annotation$genes$gene.name),]
GRN@annotation$genes[GRN@annotation$genes$gene.name == "Cd244a" & !is.na(GRN@annotation$genes$gene.name),]

# no peak anno correction needed because no GC data pulled

# Correct Gene Annotation ####

#summarized changes: force biomart with selected version AND changed host to be https://nov2020.archive.ensembl.org
populateGeneAnnotation_CUSTOM = function(GRN, EnsemblVersion = "102") {
  
  
  countsRNA.m  = getCounts(GRN, type = "rna", permuted = FALSE, asMatrix = TRUE, includeFiltered = TRUE)
  
  futile.logger::flog.info(paste0(" Calculate statistics for each of the ", nrow(countsRNA.m), " genes that were provided with the RNA-seq data (mean and CV)"))
  
  
  rowMeans_rna = rowMeans(countsRNA.m)
  rowMedians_rna = matrixStats::rowMedians(countsRNA.m)
  CV_rna = matrixStats::rowSds(countsRNA.m) /  rowMeans_rna
  
  genomeAnnotation.df = retrieveAnnotationData_CUSTOM(GRN@config$parameters$genomeAssembly, EnsemblVersion = EnsemblVersion)
  
  metadata_rna = tibble::tibble(gene.ENSEMBL = rownames(countsRNA.m), 
                                gene.mean = rowMeans_rna, 
                                gene.median = rowMedians_rna, 
                                gene.CV = CV_rna) %>%
    dplyr::left_join(genomeAnnotation.df, by = c("gene.ENSEMBL")) # %>%
  # dplyr::mutate(gene.type = forcats::fct_na_value_to_level(.data$gene.type, level = "unknown/missing"))
  
  GRN@annotation$genes = metadata_rna
  
  GRN
  
}
retrieveAnnotationData_CUSTOM = function(genomeAssembly, EnsemblVersion = NULL) {
  
  
  
  #if (source == "biomaRt") {
  
  futile.logger::flog.info(paste0("Retrieving genome annotation data from biomaRt for ", genomeAssembly, "... This may take a while"))
  
  params.l = .getBiomartParameters(genomeAssembly, suffix = "_gene_ensembl")
  
  
  columnsToRetrieve = c("chromosome_name", "start_position", "end_position",
                        "strand", "ensembl_gene_id", "gene_biotype", "external_gene_name")
  
  ensembl = .biomart_getEnsembl(biomart = "genes", version = EnsemblVersion, host = params.l[["host"]],  dataset = params.l[["dataset"]])
  
  results.df = .callBiomart(mart =  ensembl, attributes = columnsToRetrieve) 
  
  genes.df = results.df %>%
    tibble::as_tibble() %>%
    dplyr::filter(stringr::str_length(.data$chromosome_name) <= 5) %>%
    dplyr::mutate(chromosome_name = paste0("chr", .data$chromosome_name)) %>%
    dplyr::rename(gene.chr = "chromosome_name", gene.start = "start_position", gene.end = "end_position", 
                  gene.strand = "strand", gene.ENSEMBL = "ensembl_gene_id", gene.type = "gene_biotype", gene.name = "external_gene_name") %>%
    tidyr::replace_na(list(gene.type = "unknown")) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(gene.type = dplyr::recode_factor(.data$gene.type, lncRNA = "lincRNA")) %>%  # there seems to be a name change from lincRNA -> lncRNA, lets change it here 
    dplyr::mutate(gene.strand = factor(.data$gene.strand, levels = c(1,-1), labels = c("+", "-")))
  
  
  #}
  
  genes.df
  
}
.getBiomartParameters <- function(genomeAssembly, suffix = "") {
  
  host = "https://nov2020.archive.ensembl.org" # "https://www.ensembl.org"
  
  if (grepl(x = genomeAssembly, pattern = "^hg\\d+")) {
    dataset = "hsapiens"
    if (genomeAssembly == "hg38") {
    } else if (genomeAssembly == "hg19") {
      host = "https://grch37.ensembl.org"
    }
  } else if (grepl(x = genomeAssembly, pattern = "^mm\\d+")) {
    dataset = "mmusculus"
  } else if (grepl(x = genomeAssembly, pattern = "^rn\\d+")) {
    dataset = "rnorvegicus"
  } else if (grepl(x = genomeAssembly, pattern = "^dm\\d+")) {
    dataset = "dmelanogaster"
  }
  
  
  
  list(dataset = paste0(dataset, suffix), host = host)
}
.biomart_getEnsembl <- function(biomart, version, host, dataset, maxAttempts = 40) {
  
  ensembl <- NULL
  mirrorIndex <- 0
  attemptsDone = 0
  
  mirrors = c('www', 'uswest', 'useast', 'asia')
  while (!"Mart" %in% class(ensembl) && attemptsDone <= maxAttempts ) {
    mirrorIndex <- (mirrorIndex %% 4) + 1
    attemptsDone = attemptsDone + 1
    
    ensembl = tryCatch({ 
      biomaRt::useEnsembl(biomart = biomart , version = version, host = host,  
                          dataset = dataset, mirror = mirrors[mirrorIndex])
      
      
    }, error = function(e) {
      futile.logger::flog.warn(paste0("Attempt ", attemptsDone, " out of ", maxAttempts, " failed. Another automatic attempt will be performed using a different mirror. The error message was: ", e))
    })
  } 
  
  if (!"Mart" %in% class(ensembl)) {
    
    error_Biomart = paste0("A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl despite trying ", 
                           maxAttempts, 
                           " times via different mirrors. This is often caused by an unresponsive Ensembl site.", 
                           " Try again at a later time. Note that this error is not caused by GRaNIE but external services.")
   # .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = FALSE)
    return(NULL)
    
  } 
  
  futile.logger::flog.info(paste0(" Retrieving BioMart database succeeded"))
  
  ensembl
}
.callBiomart <- function(mart, attributes = NULL, values = "", filters = "", maxAttempts = 40) {
  
  result.df <- NULL
  attemptsDone = 0
  
  while (!is.data.frame(result.df) && attemptsDone <= maxAttempts ) {
    attemptsDone = attemptsDone + 1
    
    result.df = tryCatch({ 
      biomaRt::getBM(attributes = attributes,
                     filters = filters,
                     values = values,
                     mart = mart) 
      
    }, error = function(e) {
      futile.logger::flog.warn(paste0("Attempt ", attemptsDone, " out of ", maxAttempts, " failed. Another automatic attempt will be performed using a different mirror. The error message was: ", e))
    })
  } 
  
  if (!is.data.frame(result.df)) {
    
    error_Biomart = paste0("A temporary error occured with biomaRt::getBM or biomaRt::useEnsembl despite trying ", 
                           maxAttempts, 
                           " times via different mirrors. This is often caused by an unresponsive Ensembl site.", 
                           " Try again at a later time. Note that this error is not caused by GRaNIE but external services.")
   # .checkAndLogWarningsAndErrors(NULL, error_Biomart, isWarning = FALSE)
    return(NULL)
    
  } 
  
  futile.logger::flog.info(paste0(" Retrieving genome annotation succeeded"))
  
  result.df
  
  
}


GRN = populateGeneAnnotation_CUSTOM(GRN, "102")

GRN@annotation$genes[GRN@annotation$genes$gene.name == "Itgav" & !is.na(GRN@annotation$genes$gene.name),]
GRN@annotation$genes[GRN@annotation$genes$gene.name == "Icos" & !is.na(GRN@annotation$genes$gene.name),]
GRN@annotation$genes[GRN@annotation$genes$gene.name == "Cd244a" & !is.na(GRN@annotation$genes$gene.name),]


# plot QC ####
GRN = plotPCA_all(GRN, data = c("rna","peaks"), topn = 500, type = "normalized", plotAsPDF = T,
                  #pages = c(2, 3, 14),
                  forceRerun = TRUE)

# Add TF and Motif Information ####

motifs_path = "/Scottbrowne/members/smd/motifs/HOCOMOCO/v11/GRaNIE"

GRN = addTFBS(GRN,
              motifFolder = motifs_path,
              translationTable = "translationTable.tsv",
              translationTable_sep = "\t",
              TFs = "all",
              filesTFBSPattern = "",
              fileEnding = ".bed", 
              forceRerun = TRUE)

GRN = overlapPeaksAndTFBS(GRN, nCores = 1, forceRerun = TRUE)

# TF to peak connections ####

GRN = addConnections_TF_peak(GRN,
                             plotDiagnosticPlots = TRUE,
                             connectionTypes = c("expression"),
                             corMethod = "pearson",
                             forceRerun = TRUE)
# check number of connections
GRN@connections$TF_peaks[["0"]]$main %>% nrow()
# QC plots
GRN = plotCorrelations(GRN,
                       type = "TF-peak",
                       #nMax = 1,
                       #min_abs_r = 0.8, 
                       plotsPerPage = c(1,1),
                       plotAsPDF = TRUE)
GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = TRUE) #, pages = c(1)) 
GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("background"), plotAsPDF = TRUE) #, pages = c(1)) 

saveRDS(GRN, paste0(proj_path,"/data/GRANIE/", Sys.Date(),"_TF-peak_pearson_GRN.rds"))

# peak to gene connections ####

GRN = addConnections_peak_gene(GRN, 
                               corMethod = "pearson", 
                               overlapTypeGene = "TSS",
                               promoterRange = 250000, # this is the search region to make peak-gene links NOT excluding as promoter
                               TADs = NULL, 
                               #knownLinks = , # should I add the -23.8kb Pdcd1 enh?
                               nCores = 8, 
                               plotDiagnosticPlots = TRUE, 
                               plotGeneTypes = list(c("protein_coding")), # "all"
                               forceRerun = TRUE)

saveRDS(GRN, paste0(proj_path,"/data/GRANIE/",
                     Sys.Date(),"_peak-gene_pearson_GRN.rds"))

#GRN = readRDS(paste0(proj_path,"/data/GRANIE/",
#                     "2024-09-01","_peak-gene_pearson_GRN.rds")) 

# Combine Connections and Filter #####

GRN = filterGRNAndConnectGenes(GRN, 
                               TF_peak.fdr.threshold = 0.3, #0.2 # FDR not correlation R
                               # do we want to ignore TF-peak FDR given flaw of expression use?
                               # plus it's already pre-selected for FDR < 0.3 from the first steps of the workflow
                               peak_gene.fdr.threshold = 0.2, #0.2 # FDR not correlation R
                               peak_gene.fdr.method = "BH", 
                               gene.types = c("protein_coding"),
                               #peak_gene.r_range = c(0.5, 1),
                               allowMissingTFs = FALSE,
                               allowMissingGenes = FALSE,
                               forceRerun = TRUE)

# TF to gene connections #####
GRN = add_TF_gene_correlation(GRN, 
                              corMethod = "pearson", 
                              nCores = 4, forceRerun = TRUE)
saveRDS(GRN, paste0(proj_path,"/data/GRANIE/",
                    Sys.Date(),"_TF-peak-gene_pearson_GRN.rds"))

# save peak to gene ####

p2g =  getGRNConnections(GRN, 
                         type = "peak_genes", # "TF_peaks", "peak_genes", "TF_genes", "all", "all_filtered"
                         #include_TF_gene_correlations = TRUE,
                         #include_TFMetadata = TRUE,
                         include_peakMetadata = TRUE,
                         include_geneMetadata = TRUE)
write.table(p2g, file = paste0(proj_path,"/data/GRANIE/",
                                Sys.Date(),"_peak2gene_pearson.tsv"),
             quote = F,
             row.names = F,
             sep = "\t")
 
 saveRDS(p2g, file = paste0(proj_path,"/data/GRANIE/",
                            Sys.Date(),"_peak2gene_pearson.rds"))

# get eGRN df ####

GRN_connections.all = getGRNConnections(GRN, 
                                        type = "all.filtered", # "TF_peaks", "peak_genes", "TF_genes", "all", "all_filtered"
                                        include_TF_gene_correlations = TRUE,
                                        include_TFMetadata = TRUE,
                                        include_peakMetadata = TRUE,
                                        include_geneMetadata = TRUE)

GRN_connections.all

write.table(GRN_connections.all, paste0(proj_path, "/data/GRANIE/",
                                        Sys.Date(), "_GRN-connections_TFpeak0.3FDR_peakgene0.2FDR_pearson.tsv"),
            row.names = F, quote = F, sep = "\t")
saveRDS(GRN_connections.all, paste0(proj_path, "/data/GRANIE/",
                                    Sys.Date(), "_GRN-connections_TFpeak0.3FDR_peakgene0.2FDR_pearson.rds"))

GRN_connections.all[GRN_connections.all$peak_gene.r > 0.7,]$gene.name %>% unique %>% sort
GRN_connections.all[GRN_connections.all$gene.name == "Itgav" &
                      GRN_connections.all$peak_gene.r > 0.5 &
                      GRN_connections.all$TF_peak.r > 0.3,
                    c('TF.ID','peak.ID','peak_gene.distance',
                      'peak_gene.r', 'TF_peak.r', 'TF_gene.r',
                      'TF_peak.fdr','peak_gene.p_adj','TF_gene.p_raw')] # Itgav, Havcr2, Cd244a, Slamf7, Slamf1, Cd48, Ptpn22
GRN_connections.all[GRN_connections.all$gene.name == "Itgav",] # Itgav, Havcr2, Cd244a, Slamf7, Slamf1, Cd48, Ptpn22
GRN_connections.all$gene.name %>% unique() %>% sort() %>% grep(pattern = "Ptpn", value = T)

# ####

GRN = generateStatsSummary(GRN, TF_peak.fdr = c(0.05, 0.1, 0.2), TF_peak.connectionTypes = "all",
                           peak_gene.fdr = c(0.1, 0.2), peak_gene.r_range = c(0, 1),
                           allowMissingGenes = c(FALSE,
                                                 TRUE), 
                           allowMissingTFs = c(FALSE), 
                           gene.types = c("protein_coding", "lincRNA"),
                           forceRerun = TRUE)


GRN = plot_stats_connectionSummary(GRN, type = "heatmap", plotAsPDF = TRUE) # boxplot
