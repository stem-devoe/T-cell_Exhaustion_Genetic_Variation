library(dplyr)
library(stringr)
library(ggplot2)
# updated 20240903 to use counts filtered decide test from phenotype comparisons (within same strain)


source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-05_Functions_MotifDisruption_2sample.R")
motif_db = "HOCOMOCO" 
phenotype = "Naive"

save_date = "2024-09-03"
motifs = list.files(path = "/scratch2/devoes/SD029/HOCOMOCO/categorized",
                    pattern = ".bed") %>%
  str_match(pattern = "HOCOMOCO_(.*)_categorized") %>% .[,2]


# assemble df for anova ####

# load tf.df from SD029-04
tf.df = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                       motif_db, 
                       "/150bp/padj_0.1/FC_1.5/",
                       phenotype,
                       "/",
                       save_date,
                       "_tf-df.rds")
) 

# test LFC for disturbed sites ####

# run
test.spec = "B6"
# test all motifs
lfc.res = data.frame(motif = motifs,
                     p.value = lapply(motifs, 
                                      FUN = testLFC_2sample, 
                                      tf.df, 
                                      test.spec, 
                                      min_sites = 5,
                                      alt = "less") %>% 
                       unlist())
# multiple testing corrections
lfc.res$p.adj = p.adjust(lfc.res$p.value, method = "fdr")
# get means
lfc.res$mean_B6_lfc = lapply(motifs, FUN = getMeanLFC, tf.df, test.spec) %>% unlist()
lfc.res$mean_shared_lfc = lapply(motifs, FUN = getMeanLFC, tf.df, "Shared") %>% unlist()
# get difference in mean
lfc.res$diff_mean_LFC = lfc.res$mean_B6_lfc - lfc.res$mean_shared_lfc
# set labels for plotting
lfc.res$label = NA
label_bool = (lfc.res$diff_mean_LFC < -0.5 ) & (lfc.res$p.adj < 0.1)
lfc.res$label[label_bool] = lfc.res$motif[label_bool]

lfc.res %>% arrange(desc(diff_mean_LFC)) %>% tail(n = 10)
lfc.res[lfc.res$p.adj < 0.1, ]
lfc.res[grepl("NFAT",lfc.res$motif, ignore.case = T),]


saveRDS(lfc.res, 
        file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                      motif_db, 
                      "/150bp/padj_0.1/FC_1.5/",
                      phenotype,
                      "/", 
                      Sys.Date(),
                      "_lfc-res_2sample.rds")
)


# plot volcano
p.volcano = ggplot(data = lfc.res, aes(y= -log10(p.adj), x = diff_mean_LFC, label = label)) + 
  geom_point() + 
  ggrepel::geom_text_repel(min.segment.length = 0, color = "purple") +
  geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), color = "red", linetype = "dashed") +
  theme_classic() +
  xlab(paste0("difference in mean LFC"))

pdf(paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
           motif_db, 
           "/150bp/padj_0.1/FC_1.5/",
           phenotype,
           "/",
           Sys.Date(),
           "_2sample_MotifDisruption_volcano.pdf"))
p.volcano
dev.off()


# plot rank #####

p.waterfall = lfc.res %>% 
  arrange(desc(diff_mean_LFC)) %>% tibble::add_column(rank = seq_len(nrow(.))) %>%
  ggplot(aes(y = diff_mean_LFC, x = rank, label = motif)) + geom_point() + 
  ggrepel::geom_text_repel(data = ~tail(.x, 10), max.overlaps = 20, min.segment.length = 0, color = "purple") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "red") +
  theme_classic()

pdf(file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                  motif_db,
                  "/150bp/padj_0.1/FC_1.5/",
                  phenotype,
                  "/",
                  Sys.Date(),
                  "_2sample_MotifDisruptions_LFC_waterfall.pdf"))
p.waterfall
dev.off()

# plot scatter ####

p.scatter = lfc.res %>%
  ggplot(aes(y = mean_B6_lfc, x = mean_shared_lfc)) +
  geom_point(aes(color = p.adj), size = 3, shape = 16) +  scale_color_viridis_c(direction = -1) +
  ggrepel::geom_text_repel(aes(label = label),min.segment.length = 0, color = "black") +
  theme_minimal()

pdf(file = paste0("/Scottbrowne/members/smd/Projects/SD029/motif_disruption/",
                  motif_db,
                  "/150bp/padj_0.1/FC_1.5/",
                  phenotype,
                  "/",
                  Sys.Date(),
                  "_2sample_MotifDisruptions_LFC_scatter.pdf"))
p.scatter
dev.off()
