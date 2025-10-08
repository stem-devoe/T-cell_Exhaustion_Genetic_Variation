install.packages("devtools")
devtools::install_gitlab("grp-zaugg/GRaNIE", host = "git.embl.de", subdir = "src/GRaNIE", force = TRUE)
BiocManager::install("Rsubread")
BiocManager::install(c("org.Mm.eg.db", "BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
remotes::install_version("dbplyr", version = "2.3.4" ) # https://github.com/Bioconductor/AnnotationHub/issues/46
BiocManager::install("edgeR")
