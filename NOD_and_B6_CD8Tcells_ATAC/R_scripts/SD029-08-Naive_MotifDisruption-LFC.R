library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(limma)
# updated 20240903 to use counts filtered decide test from phenotype comparisons (within same strain)


set.seed(901)

source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-04-Functions_MotifDisruption-LFC.R")

# load metadata ####

phenotype = "Naive" 
motif_db = "HOCOMOCO" 

if(!dir.exists(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                      motif_db,
                      "/150bp/padj_0.1/FC_1.5/",
                      phenotype))){
  dir.create(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                    motif_db,
                    "/150bp/padj_0.1/FC_1.5/",
                    phenotype), recursive = T)
}

meta = read.delim("/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv")


# load limma-voom tfit for strain comparisons ####

tfit.strain = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/150bp/tfit.rds")
tt.strain.naive = topTreat(tfit.strain, coef = "NOD-B6_Naive", n = Inf)


# 1. Define peakset ####

# load consensus peaks 
peaks = rtracklayer::import("/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.renamed.bed")


# categorize phenoype comparisons
# load counts filtered decide tests results of phenotype comparisons
pheno.res = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/decideTests-Phenotype_FC1.5_countsFiltered.rds")
cate.pheno.cart = apply(pheno.res[,c('CART-Naive.NOD', 'CART-Naive.B6')], 
                        1, 
                        FUN = categorizePeakChange)
cate.pheno.lcmv = apply(pheno.res[,c('LCMV-Naive.NOD', 'LCMV-Naive.B6')], 
                        1, 
                        FUN = categorizePeakChange)
table(cate.pheno.cart)
table(cate.pheno.lcmv)

cate.pheno.cart = data.frame(peakID = names(cate.pheno.cart), category = cate.pheno.cart)
cate.pheno.lcmv = data.frame(peakID = names(cate.pheno.lcmv), category = cate.pheno.lcmv)

cate.pheno = inner_join(cate.pheno.cart, cate.pheno.lcmv, 
                        by = "peakID", 
                        suffix = c("CART","LCMV") )

# subset peaks for peaks LOST in both CART AND LCMV 

#cate.pheno[grepl("L", cate.pheno$categoryCART) & grepl("L",cate.pheno$categoryLCMV),]
peak_subset_gr = peaks[peaks$name %in% cate.pheno$peakID[grepl("L", cate.pheno$categoryCART) & grepl("L",cate.pheno$categoryLCMV)]] # 

saveRDS(peak_subset_gr,
        paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
               motif_db,
               "/150bp/padj_0.1/FC_1.5/",
               phenotype,
               "/Lost_CART-AND-LCMV_to-Naive_any-strain_peak_subset.rds"))
rtracklayer::export(peak_subset_gr,
                    paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                           motif_db,
                           "/150bp/padj_0.1/FC_1.5/",
                           phenotype, 
                           "/Lost_CART-AND-LCMV_to-Naive_any-strain_peak_subset.bed"))
# 2. categorize motif sites within peaks ####

# run 
motifs = list.files(path = paste0("/scratch2/devoes/SD029/",motif_db,"/categorized"),
                    pattern = ".bed") %>%
  str_match(pattern = paste0(motif_db,"_(.*)_categorized")) %>% .[,2]

cate.motif = getMotifCategories(motif_names = motifs,
                                motif_db = motif_db,
                                motif_path = paste0("/scratch2/devoes/SD029/",
                                                    motif_db,
                                                    "/categorized"),
                                peakset_gr = peak_subset_gr)


saveRDS(cate.motif, paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                           motif_db,
                           "/150bp/padj_0.1/FC_1.5/",
                           phenotype,
                           "/",
                           Sys.Date(),
                           "_LFC_cate-motif_",
                           motif_db,
                           ".rds"))
# cate.motif = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
#                             motif_db,
#                             "/150bp/padj_0.1/FC_1.5/",
#                             phenotype,
#                             "/",
#                             save_date, # Sys.Date part of file # 2024-07-22
#                             "_LFC_cate-motif_",motif_db,".rds"))

# assemble df for anova ####

tf.df = left_join(cate.motif, tt.strain.naive, by = "GeneID") # dont need to pre-order to match if using by in left join

# test LFC for disturbed sites ####

# run
test.spec = "B6"
lfc.res = data.frame(motif = motifs,
                     p.value = lapply(motifs, 
                                      FUN = testLFC, 
                                      tf.df, 
                                      test.spec, 
                                      min_sites = 5,
                                      alt = "less") %>% 
                       unlist())
# multiple testing corrections
lfc.res$p.adj = p.adjust(lfc.res$p.value, method = "fdr")
# get means
lfc.res$mean_lfc = lapply(motifs, FUN = getMeanLFC, tf.df, test.spec) %>% unlist()
# set labels for plotting
lfc.res$label = NA
label_bool = (abs(lfc.res$mean_lfc) > log2(1.5)) & (lfc.res$p.adj < 0.1)
lfc.res$label[label_bool] = lfc.res$motif[label_bool]

lfc.res %>% arrange(desc(mean_lfc)) %>% tail(n = 10)
lfc.res[lfc.res$p.adj < 0.1, ]

saveRDS(lfc.res, 
        file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                      motif_db,
                      "/150bp/padj_0.1/FC_1.5/",
                      phenotype,
                      "/", 
                      Sys.Date(),
                      "_lfc-res.rds")
)
saveRDS(tf.df,
        file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                      motif_db,
                      "/150bp/padj_0.1/FC_1.5/",
                      phenotype,
                      "/",
                      Sys.Date(),
                      "_tf-df.rds")
)

# plot volcano
p.volcano = ggplot(data = lfc.res, aes(y= -log10(p.adj), x = mean_lfc, label = label)) + 
  geom_point() + 
  ggrepel::geom_text_repel(min.segment.length = 0, color = "purple") +
  geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), color = "red", linetype = "dashed") +
  theme_classic() +
  xlab(paste0("mean L2FC for ", test.spec, " Specific Sites"))

pdf(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
           motif_db,
           "/150bp/padj_0.1/FC_1.5/",
           phenotype,
           "/",
           Sys.Date(),
           "_MotifDisruption_volcano.pdf"))
p.volcano
dev.off()


# plot rank #####

p.waterfall = lfc.res %>% 
  arrange(desc(mean_lfc)) %>% tibble::add_column(rank = seq_len(nrow(.))) %>%
  ggplot(aes(y = mean_lfc, x = rank, label = motif)) + geom_point() + 
  ggrepel::geom_text_repel(data = ~tail(.x, 10), max.overlaps = 20, min.segment.length = 0, color = "purple") +
  geom_hline(yintercept = -log2(1.5), linetype = "dashed", color = "red") +
  theme_classic()

pdf(file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                  motif_db,
                  "/150bp/padj_0.1/FC_1.5/",
                  phenotype,
                  "/",
                  Sys.Date(),
                  "_MotifDisruptions_LFC_waterfall.pdf"))
p.waterfall
dev.off()

# plot distributions ####

p.tfa = lapply(motifs, FUN = plotTFA, tf.df, lfc.res)
p.tfa[[2]]

pdf(file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                  motif_db,
                  "/150bp/padj_0.1/FC_1.5/",
                  phenotype,
                  "/",
                  Sys.Date(),
                  "_MotifDisruptions_LFC.pdf"))
p.tfa
dev.off()
