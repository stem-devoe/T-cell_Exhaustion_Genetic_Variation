library(ggplot2)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(Rsubread)
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

fold_change = c(1.5,2,3)
padj = 0.1

group_colors =  c(B6_Naive = "#ff33ff",
                  NOD_Naive = "#cc99ff",
                  B6_Memory = "#ff6600",
                  NOD_Memory = "#ff9999",
                  B6_LCMVcl13 = "#009900",
                  NOD_LCMVcl13 = "#99cc66",
                  #B6_CART = "#3333ff",
                  B6_CART = "#4969ff",
                  NOD_CART = "#6699ff")


# 3. load data ####

# load sample data

metadata = read.delim(metadata_file)

# set bam paths (NOD and B6 are different)

bam_files = paste0(data_path,
                   "/",
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

# see /Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-01_lmv_countsFilter.R

counts = readRDS(paste0(out_data, "/ArchR",half_width,"bpExt_featureCounts.rds"))

# 5. filter peaks by cpm ####

# get library size 
lib_size = colSums(counts$counts)
quantile(lib_size, probs = seq(0,1,0.1)) 

# get cpm per sample per peak
counts_per_million = cpm(counts$counts)
colnames(counts_per_million) = colnames(counts_per_million) %>% str_remove(pattern = "_REP1.*")
if(!all(metadata$SampleID == colnames(counts_per_million))){
  counts_per_million = counts_per_million[,metadata$SampleID]
} else{
  print("order matches metadata")
}

# get avg cpm by group and overall avg cpm
cpm_group = data.frame(t(counts_per_million), group = metadata$Group)
cpm_group = cpm_group %>% group_by(group) %>% summarise(across(.cols = everything(), mean))

cpm_avg = colMeans(cpm_group[,!grepl('group',colnames(cpm_group))])

# check cpm distributions
quantile(cpm_group[,!grepl('group',colnames(cpm_group))] %>% as.vector() %>% unlist(), 
         probs = seq(0,1,0.05))

quantile(cpm_avg, probs = seq(0,1,0.1))
quantile(cpm_avg, probs = seq(0,1,0.25))
hist(cpm_avg[which(cpm_avg < 30)], breaks = 200) 
abline(v=quantile(cpm_avg, .1), col='red', lwd=2)
abline(v=quantile(cpm_avg, .05), col='orange', lwd=2)
abline(v=1, col='blue', lwd=2)

# set cpm_filter
cpm_filt = 3 # a little more lenient than 5th percentile, smallest lib would be ~7 counts, largest ~45
num_groups_cpm_filt = 1

# want peak to be prominent in at least one group
keep_peaks = colSums(cpm_group[,!grepl('group',colnames(cpm_group))] > cpm_filt) >= num_groups_cpm_filt
table(keep_peaks)
keep_peaks = names(keep_peaks)[keep_peaks]
head(keep_peaks)

# filter
counts_pass = counts$counts[keep_peaks,]
head(peak_saf)
rownames(peak_saf) = peak_saf$GeneID
peak_pass = peak_saf[keep_peaks,]

# 6. build design ####

# set up DGEList object to store counts matrix and metadata variables for model
# this is from EdgeR but limma-voom is built on it
y = DGEList(counts = counts_pass,
            group = metadata$Group,
            genes = peak_pass,
            samples = metadata) 

# set reference state
y$samples$Background <- relevel(as.factor(y$samples$Background), ref = "B6")
y$samples$Phenotype <- relevel(as.factor(y$samples$Phenotype), ref = "Naive")

# voom
design <- model.matrix(~Background*Phenotype, y$samples) # each factor is a variable, not a covariate, so with and without intercept should be equivalent models
print(design)

#v <- voom(y, design, plot=TRUE, normalize="quantile")

v = readRDS( paste0(out_data,"/v_",half_width,"bp.rds"))

# 8. run limma model ####

fit <- lmFit(v, design)

# define contrasts for comparison: NOD - B6
c_CART_B6 = c(0,0,1,0,0,0) 
c_LCMV_B6 = c(0,0,0,1,0,0) 
c_CART_NOD = c(0,0,1,0,1,0) 
c_LCMV_NOD = c(0,0,0,1,0,1) 

contrastMX <- matrix(c(c_CART_B6,c_LCMV_B6,c_CART_NOD,c_LCMV_NOD),nrow=6)
colnames(contrastMX) <- c("CART-Naive.B6","LCMV-Naive.B6","CART-Naive.NOD","LCMV-Naive.NOD")

# compute estimate coefficients and standard error for each contrast from the linear model
cfit <- contrasts.fit(fit, contrasts = contrastMX)


tfit <- list()
res = list()
for(i in 1:length(fold_change)){
  tfit[[i]] <- treat(cfit, lfc = 0) 
  # use decideTests for inclusion of post-hoc. topTreat does an ANOVA F-test which does not test a specific contrast rather if any coefficient is non-zero
  ## NOTE : decideTests REORDERS the peak regions! ##
  res[[i]] <- decideTests(tfit[[i]], 
                          adjust.method = "BH", 
                          lfc = log2(fold_change[i]),
                          p.value = padj)
  # sort results to match peak order of filtered saf
  res[[i]] = res[[i]][match(peak_pass$GeneID, rownames(res[[i]])),]
  print(colSums(res[[i]] != 0))
  print(length(which(rowSums(res[[i]] != 0) != 0)))
}

# check correlation of L2FC to pvalues for validity of filtering by L2FC in decideTests
df.tt = lapply(seq_len(ncol(contrastMX)), FUN = topTreat, fit = tfit[[2]], number = 1e6) # CART coef = 1, lcmv = 2 , mem = 3, naive = 4
p.volc = list()
for(i in seq_len(ncol(contrastMX))){p.volc[[i]] = ggplot(data = df.tt[[i]], aes(x = logFC, y = -log10(P.Value))) + ArchR::theme_ArchR() + geom_point(size = 1, shape = 16)  + ggtitle(colnames(tfit[[2]]$coefficients)[i])}
gridExtra::grid.arrange(grobs = p.volc, ncol = 2)
for(i in seq_len(ncol(contrastMX))){
  print(colnames(tfit[[2]]$coefficients)[i])
  print(paste0("Abs L2FC to P Value correlation: ", cor(x = abs(df.tt[[i]]$logFC), y = df.tt[[i]]$P.Value)))
  print(paste0("Abs L2FC to -log10Pval correlation: ", cor(x = abs(df.tt[[i]]$logFC), y = -log10(df.tt[[i]]$P.Value))))
}

# write decideTests to file
for(i in 1:length(fold_change)){
  
  if(!dir.exists(paste0(out_data,"/",half_width,"bp/padj_",padj))){
    print("create output directories")
    dir.create(paste0(out_data,"/",half_width,"bp/padj_",padj), recursive = T)
  }
  
  saveRDS(res[[i]], 
          file = paste0(out_data,"/",half_width,"bp/padj_",padj,"/decideTests-Phenotype_FC",fold_change[i],".rds"))
}

# write tfit to file - even though I had the loop, tfit is the same at this stage it is generated

saveRDS(tfit[[i]], 
        file = paste0(out_data,"/",half_width,"bp/tfit_phenotype.rds"))

# 9. Filter for minimum signal ####

min_signal_counts = 25 # approx 5 cpm 

# get average signal per group
# if i come back to clean this up switch to group_by and summarize(across())
counts_list = split.data.frame(t(counts_pass), f = metadata$Group) %>% 
  lapply(., FUN = t)
avg_cpm_pass = lapply(counts_list, FUN = rowMeans) %>% rlist::list.cbind() %>% as.data.frame()


# check that peak orders are the same for counts mtx and results from decideTests
lapply(res, FUN = function(x, counts.mtx){all(rownames(x) == rownames(counts.mtx))}, counts.mtx = counts_pass)

# for each phenotype 
compari = colnames(contrastMX)
res.countsFilt = list()
for(j in 1:length(fold_change)){
  temp = {}
  for(i in seq_len(ncol(contrastMX))){
    #
    geno = str_remove(colnames(contrastMX)[i], pattern = ".*-Naive.")
    pheno = str_remove(colnames(contrastMX)[i], pattern = "-Naive.*")
    search_term = c(paste0(geno,"_",pheno), paste0(geno,"_Naive"))
    search_term = paste0(search_term, collapse = "|")
    # get peaks passing count filter per group
    count_bool = rowSums(avg_cpm_pass[,grep(search_term, colnames(avg_cpm_pass))] > min_signal_counts) > 0 # should we check for min count in Exh Pehno and Naive or JUST Exh Pheno
    ## determine filter passing and significantly gained in NOD or B6
    # get ind of sig peak from decideTests for counts mtx
    sig_up_bool = rownames(counts_pass) %in% rownames(res[[j]][,colnames(contrastMX)[i]])[which(res[[j]][,colnames(contrastMX)[i]] == 1)]
    sig_down_bool = rownames(counts_pass) %in% rownames(res[[j]][,colnames(contrastMX)[i]])[which(res[[j]][,colnames(contrastMX)[i]] == -1)]
    # determine which peaks pass signal filter AND were sig from decide Tests
    up_ind = which(count_bool & sig_up_bool)
    down_ind = which(count_bool & sig_down_bool)
    sig_ind = rep(0, nrow(counts_pass))
    sig_ind[up_ind] = 1
    sig_ind[down_ind] = -1
    temp = cbind(temp, sig_ind)
  }
  res.countsFilt[[j]] = temp
  colnames(res.countsFilt[[j]]) = colnames(contrastMX)
  rownames(res.countsFilt[[j]]) = rownames(counts_pass)
  print(colSums(res.countsFilt[[j]] != 0))
  print(length(which(rowSums(res.countsFilt[[j]] != 0) != 0)))
  
}


# write decideTests to file
for(i in 1:length(fold_change)){
  saveRDS(res.countsFilt[[i]], 
          file = paste0(out_data,"/",half_width,"bp/padj_",padj,"/decideTests-Phenotype_FC",fold_change[i],"_countsFiltered.rds"))
}

# 10. write results ####


# write significant, signal passing peaks to bed files
# all, gained in nod, gained in b6

writeSigPeaks <- function(results, peaks_in_model, out_base, write_all = T){
  
  # convert saf to bed
  peak_bed = peaks_in_model[,c(2,3,4,1)]
  peak_bed[,2] = peak_bed[,2] - 1 # convert from 1 to 0 based
  
  # sort results to match order of peak_saf
  res_sort_ind = match(peaks_in_model$GeneID, rownames(results))
  results = results[res_sort_ind,]
  
  for(i in 1:ncol(results)){
    
    # get name of comparison
    contrast = colnames(results)[i]
    
    # set file path name to save results
    out_bed = paste0(out_base,"/",contrast)
    
    # get index of significant peaks
    ind_up = which(results[,i] == 1)
    ind_down = which(results[,i] == -1)
    
    # get significant peaks
    peak_up = peak_bed[ind_up,]
    peak_down = peak_bed[ind_down,]
    
    # write peaks to bed file
    write.table(peak_up,
                file = paste0(out_bed,".up.bed"),
                sep = "\t",
                col.names = F,
                row.names = F,
                quote = F)
    
    write.table(peak_down,
                file = paste0(out_bed,".down.bed"),
                sep = "\t",
                col.names = F,
                row.names = F,
                quote = F)
    
    # save all significant peaks regardless of where accessibility is gained or lost
    if(write_all){
      
      ind_sig = sort(unique(c(ind_up,ind_down)))
      peak_all = peak_bed[ind_sig,]
      
      write.table(peak_all,
                  file = paste0(out_bed,".all.bed"),
                  sep = "\t",
                  col.names = F,
                  row.names = F,
                  quote = F)
      
    }
    
  }
  

  
  return()
  
}

# check that order of peaks in voom model is the same as cpm filter passing count matrix
all(rownames(counts_pass) == v$genes$GeneID)
all(rownames(res.countsFilt[[1]]) == v$genes$GeneID)

for(i in 1:length(fold_change)){
  if(!dir.exists(paste0(out_data,"/",half_width,"bp/padj_",padj,"/sig_peaks/FC_",fold_change[i]))){
    print("create output directories")
    dir.create(paste0(out_data,"/",half_width,"bp/padj_",padj,"/sig_peaks/FC_",fold_change[i]), recursive = T)
  }
  out_base = paste0(out_data,"/",half_width,"bp/padj_",padj,"/sig_peaks/FC_",fold_change[i])
  #print(out_base)
  writeSigPeaks(res.countsFilt[[i]], # save count signal passing sig peaks
                v$genes,
                out_base)
}

