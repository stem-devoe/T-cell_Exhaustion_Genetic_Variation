library(ggplot2)
library(dplyr)
library(stringr)
library(limma)
library(egg)

proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
fold_change = 1.5
half_width = 150
padj = 0.1
#pheno = "CART"

fill.colors = c(B6_Naive = "#ff33ff",
                NOD_Naive = "#cc99ff",
                B6_LCMVcl13 = "#009900",
                NOD_LCMVcl13 = "#99cc66",
                B6_CART = "#4969ff",
                NOD_CART = "#6699ff")

# scatter theme
source("/Scottbrowne/members/smd/Projects/SD029/R_scripts/SD029-PTheme.R")

dt.res = readRDS(paste0(proj.dir, "/data/",half_width, "bp/padj_",padj,"/decideTests-Phenotype_FC",fold_change,"_countsFiltered.rds"))


dt.res.summary = limma::summary.TestResults(dt.res) # i need a different summary. not base but can't remember which
#dt.res.summary = dt.res.summary[,grepl(paste0(pheno,collapse = "|"), colnames(dt.res.summary))]
dt.res.summary = reshape2::melt(dt.res.summary)
dt.res.summary = dt.res.summary[!grepl("NotSig",dt.res.summary$Var1),]
dt.res.summary$pheno = str_match(dt.res.summary$Var2, "(.*)-Naive")[,2] %>% str_replace("LCMV","LCMVcl13")
dt.res.summary$strain = str_match(dt.res.summary$Var2, "Naive.(.*)")[,2]
#dt.res.summary$Var2 = str_remove(dt.res.summary$Var2, paste0(pheno,"-Naive."))
dt.res.summary$change = str_replace(dt.res.summary$Var1, pattern = "Down",
                                  replacement = " < Naive") %>%
  str_replace(pattern = "Up", replacement = " > Naive")
dt.res.summary$change = paste0(dt.res.summary$pheno, dt.res.summary$change)
dt.res.summary$fill = fill.colors[paste0(dt.res.summary$strain,"_",dt.res.summary$pheno)]
dt.res.summary
p.bar = ggplot(dt.res.summary, aes(x = strain, y = value, fill = fill, alpha = Var1))  + 
  geom_bar(position="stack", stat="identity", width = 0.75) + 
  facet_wrap(~pheno) +
  theme_classic() + my.theme +  xaxis.gap.fix +
  theme(panel.spacing.x=unit(0, "lines")) +
  scale_fill_identity() +
  scale_alpha_manual(values = c(0.7,1)) + 
 # geom_bar_pattern(aes(pattern = Var1),width = 0.75, position = "stack",stat = "identity") + 
 # scale_pattern_manual(values = c('circle', 'none')) +
  ylab("number of peaks") +
  theme(legend.position = "none") + xlab("")
p.bar

p.bar.egg = set_panel_size(p.bar, width = unit(1, "in"), height = unit(2,"in")) # width will set for each facet, not whole plot

ggsave(paste0("/scratch2/devoes/SD029/figures/",Sys.Date(),"_bar_phenotype.pdf"),
       plot = p.bar.egg,
       device = "pdf",
       unit = "in",
       width = 3,
       height = 3,
       dpi = 300)

# 

shared.cart = vennCounts(dt.res[,grepl("CART", colnames(dt.res))], include = "up")
shared.cart.df = data.frame(change = c("B6","NOD","Shared"), 
                       num_peaks = shared.cart[2:4,3])

shared.lcmv = vennCounts(dt.res[,grepl("LCMV", colnames(dt.res))], include = "up")
shared.lcmv.df = data.frame(change = c("B6","NOD","Shared"), 
                            num_peaks = shared.lcmv[2:4,3])
