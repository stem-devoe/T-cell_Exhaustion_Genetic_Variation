library(dplyr)
library(stringr)

motif_db = "HOCOMOCO"
options(scipen = 99999)

motifs = list.files(path = paste0("/scratch2/devoes/SD029/",motif_db,"/categorized"),
                    pattern = ".bed",
                    full.names = T)
motif_names = str_match(motifs, pattern = paste0(motif_db,"_(.*)_categorized")) %>% .[,2]
names(motifs) = motif_names
motifs = lapply(motifs, FUN = function(x){y= read.delim(x, header = F, sep = "\t"); y[grepl("NOD",y[,7]),]})
motifs = lapply(motif_names, FUN = function(x, motifs){motifs[[x]]$motif = x; motifs[[x]]}, motifs)
nod = purrr::reduce(motifs, rbind)
write.table(nod, file = "/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_consensus_nodsites.bed",
            col.names = F,
            sep = "\t",
            row.names = F,
            quote = F)

