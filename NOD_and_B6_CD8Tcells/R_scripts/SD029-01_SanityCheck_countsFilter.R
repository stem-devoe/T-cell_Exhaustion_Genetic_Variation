dt = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/decideTests_FC1.5.rds")
summary(dt)
dt.filt = readRDS("/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/decideTests_FC1.5_countsFiltered.rds")
apply(dt.filt, MARGIN = 2, table)

out_data = "/Scottbrowne/members/smd/Projects/SD029/data"
half_width = 150
metadata_file = "/Scottbrowne/members/smd/Projects/SD029/metadata/low-noise_metadata.tsv"
metadata = read.delim(metadata_file)
counts = readRDS(paste0(out_data, "/ArchR",half_width,"bpExt_featureCounts.rds"))
counts_per_million = edgeR::cpm(counts$counts)
cpm_group = data.frame(t(counts_per_million), group = metadata$Group)
cpm_group = cpm_group %>% group_by(group) %>% summarise(across(.cols = everything(), mean))
#v = readRDS(paste0("/Scottbrowne/members/smd/Projects/SD029/data/v_",half_width,"bp.rds"))

cpm_lim = 5
# CART
peakIDs = rownames(dt[dt[,1] != 0,])
table(cpm_group[1,peakIDs] > cpm_lim | cpm_group[4,peakIDs] > cpm_lim) 
#LCMV
peakIDs = rownames(dt[dt[,2] != 0,])
table(cpm_group[2,peakIDs] > cpm_lim | cpm_group[5,peakIDs] > cpm_lim) 
# Naive
peakIDs = rownames(dt[dt[,3] != 0,])
table(cpm_group[3,peakIDs] > cpm_lim | cpm_group[6,peakIDs] > cpm_lim) 

apply(dt.filt, MARGIN = 2, table) %>% apply(MARGIN = 2, FUN = function(x){x[1]+x[3]})
