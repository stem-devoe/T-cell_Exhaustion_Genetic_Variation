library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(egg) # set_panel_size
#library(ggpattern)

# load data

fold_change = 1.5
half_width = 150
padj = 0.1
phenotypes = c("CART","LCMVcl13","Naive")

proj.dir = "/Scottbrowne/members/smd/Projects/SD029"
scratch.dir = "/scratch2/devoes/SD029"

# metadata
meta = read.delim(paste0(proj.dir, "/metadata/low-noise_metadata.tsv"))

# peaks
peak.consensus.file = paste0(scratch.dir,
                             "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt",
                             "/globallambda/", half_width, "bp",
                             "/All_Samples.fwp.filter.non_overlapping.bed")
peak.consensus = readRDS(paste0(scratch.dir,
                                "/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt",
                                "/globallambda/",half_width,"bp",
                                "/All_Samples.fwp.filter.non_overlapping.rds"))
peak.sig.files = list.files(paste0(proj.dir, "/data/",
                                   half_width, "bp/padj_",padj,"/sig_peaks/FC_", fold_change),
                            pattern = "up.bed|down.bed",
                            full.names = T) %>% 
 # grep(pattern = paste0(phenotypes,collapse = "|"), value = T) %>%
  grep(pattern = "NOD-B6_", value = T)
peak.sig = lapply(peak.sig.files, FUN = rtracklayer::import)
names(peak.sig) = basename(peak.sig.files) %>% str_remove(pattern = ".bed")

# counts
fcounts = readRDS(paste0(proj.dir, "/data/ArchR",half_width,"bpExt_featureCounts.rds"))


# limma voom model
lmv.v = readRDS(paste0(proj.dir, "/data/v_",half_width,"bp.rds"))

# color palettes 
sort.key = c("B6_Naive",
             "NOD_Naive",
             "B6_Memory",
             "NOD_Memory",
             "B6_LCMVcl13",
             "NOD_LCMVcl13",
             "B6_CART",
             "NOD_CART")

fill.colors = c(B6_Naive = "#ff33ff",
                NOD_Naive = "#cc99ff",
                B6_Memory = "#ff6600",
                NOD_Memory = "#ff9999",
                B6_LCMVcl13 = "#009900",
                NOD_LCMVcl13 = "#99cc66",
                #B6_CART = "#3333ff",
                B6_CART = "#4969ff",
                NOD_CART = "#6699ff")


# scatters of avg signal (counts per million) per peak between strains for each condition
pheno = unique(meta$Phenotype) #%>% grep(pattern = paste0(phenotypes,collapse = "|"), value = T)

## average peak counts by group
consensus.mtx = fcounts$counts
colnames(consensus.mtx) = colnames(consensus.mtx) %>% str_remove(pattern = "_REP1.*")
# match meta data order
#consensus.mtx = consensus.mtx[rownames(lmv.v$E), meta$SampleID] # keep only cpm passing peaks that made it to the model
consensus.mtx = consensus.mtx[, meta$SampleID] # keep all peaks

# convert to cpm to account for lib size
consensus.mtx = edgeR::cpm(consensus.mtx)
consensus.mtx.avg = aggregate(x = t(consensus.mtx),
                              by = list(meta$Group), 
                              FUN = mean) # avg cpm per group
rownames(consensus.mtx.avg) = consensus.mtx.avg[['Group.1']]
consensus.mtx.avg = consensus.mtx.avg[, -1] # remove grouping column

# get index for sig peaks  ####
sig.ind = lapply(peak.sig,
                 FUN = function(x, peak.id) {match(x$name, peak.id)},
                 peak.id = rownames(consensus.mtx))
# FUN
setColorFillsScatter = function(phenotype, cpm.mtx, sig.ind.list, sort.key, color.key){
  # isolate sig peak indices for only current phenotype of iteration
  pheno.inds = grep(names(sig.ind.list), pattern = phenotype)
  temp.l = sig.ind.list[pheno.inds]
  set.color = rep(NA, nrow(cpm.mtx))
  
  for(i in 1:length(temp.l)){
    if(grepl(names(temp.l)[i], pattern = "up")){
      set.gain = "NOD"
    } else {
      set.gain = "B6"
    }
    
    color.ind = which(grepl(sort.key, pattern = phenotype) & grepl(sort.key, pattern =  set.gain))
    set.color[temp.l[[i]]] = fill.colors[color.ind]
  }
  
  # 25 counts was used when I did the counts filtering before kmeans - 
  # given avg lib size approx 5mil: 5cpm
  
  return(set.color)
}
# RUN
scatter.fills = lapply(pheno, 
                       FUN = setColorFillsScatter,
                       cpm.mtx = consensus.mtx,
                       sig.ind.list = sig.ind,
                       sort.key = sort.key,
                       color.key = fill.colors)
names(scatter.fills) = pheno

# ID "common" peaks ####

# FUN
setCommonPeaks = function(phenotype, avg.cpm.mtx, lmv.v, FCrange = c(0.85,1.15), min.cpm = 5){
  # get condition values
  avg.cpm.mtx = t(avg.cpm.mtx)
  pheno.ind = grep(phenotype, colnames(avg.cpm.mtx))
  avg.cpm.mtx = avg.cpm.mtx[, pheno.ind]
  # keep peaks with cpm passing 
  cpm.keep = rowSums(avg.cpm.mtx >= min.cpm) == 2
  avg.cpm.mtx = avg.cpm.mtx[cpm.keep,]
  # get FC
  fc = avg.cpm.mtx[,1] / avg.cpm.mtx[,2]
  fc.keep = which(fc > FCrange[1] & fc < FCrange[2])
  # keep peaks with minimal FC
  avg.cpm.mtx = avg.cpm.mtx[fc.keep,]
  
  common.peaks = rownames(avg.cpm.mtx)
  
  return(common.peaks)
  
}

# RUN
common.peaks = lapply(pheno, 
                      FUN = setCommonPeaks,
                      avg.cpm.mtx = consensus.mtx.avg, 
                      lmv.v)
names(common.peaks) = pheno


# Plot Scatter With Common ####
# FUN
perPhenoPlotScatterWithCommon = function(phenotype, avg.cpm.mtx, common.peaks, scatter.fills, 
                                         common.color = "grey20", 
                                         common.bool = T,
                                         logtf = F, 
                                         keep.legend = F, 
                                         pt.size = 0.5, 
                                         axis.size = 7,
                                         label.size = 6,
                                         title.size = 7,
                                         border.width = 0.2){
  
  if(logtf){
    avg.cpm.mtx = t(avg.cpm.mtx) + 1 # add to avoid log(10) issues with scaling
    avg.cpm.mtx = log10(avg.cpm.mtx)
  } else{
    avg.cpm.mtx = t(avg.cpm.mtx)
  }
  
  # form df
  x.ind = which(grepl(colnames(avg.cpm.mtx), pattern = "B6") & grepl(colnames(avg.cpm.mtx), pattern = phenotype))
  #print(x.ind)
  y.ind = which(grepl(colnames(avg.cpm.mtx), pattern = "NOD") & grepl(colnames(avg.cpm.mtx), pattern = phenotype))
  #print(y.ind)
  fill.ind = grep(names(scatter.fills), pattern = phenotype)
  df = data.frame(B6 = avg.cpm.mtx[,x.ind], 
                  NOD = avg.cpm.mtx[,y.ind], 
                  fill = scatter.fills[[fill.ind]])
  
  # set color fill for common peaks
  if(common.bool){
    common.peaks = common.peaks[[grep(phenotype, names(common.peaks))]]
    common.ind = match(common.peaks, rownames(avg.cpm.mtx))
    df$fill[common.ind] = common.color
  }
  
  
  cols = factor(df$fill) %>% levels()
  
  # initiate ggplot
  p = ggplot(df, aes(x = B6, y = NOD)) + theme_classic()
  p = p + geom_point(mapping = aes(color = fill), size = pt.size, stroke = 0, shape = 16) + 
    scale_color_manual(values = cols) # set fill color
  p = p + 
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = border.width), # add border to whole plot
          axis.line = element_line(linewidth = border.width)) # set axis to same width as border
  p = p + xlab(paste0("B6 ", phenotype)) + 
    ylab(paste0("NOD ", phenotype))
  
  
  
  # format ggplot for log scaling
  if(logtf){
    p = p + ggtitle("ATAC-seq log10(cpm + 1)/peak")
  } else {
    p = p + scale_x_log10(breaks = breaks_log(), # log scale axes
                          labels = trans_format("log10", math_format(10^.x))) + 
      scale_y_log10(breaks = breaks_log(),
                    labels = trans_format("log10", math_format(10^.x))) 
    p = p + annotation_logticks(short = unit(0.025, "cm"),
                                mid = unit(0.05, "cm"),
                                long = unit(0.075, "cm"),
                                size = 0.25) # add log scale ticks
    p = p + ggtitle("ATAC-seq cpm/peak") 
  }
  
  # keep or remove legend
  if(!keep.legend){
    p = p + guides(color = "none") # remove legend  
  } else{
    p = p + theme(legend.position = "bottom", legend.title=element_blank())
  }
  
  # adjust font sizes and tick sizes
  p = p + theme(axis.title = element_text(size = axis.size, margin = margin(0,0,0,0)), # axis title
                axis.text = element_text(size = label.size, margin = margin(0,0,0,0)), # tickmark labels
                plot.title = element_text(size = title.size, margin = margin(0,0,0,0)), # plot title
                axis.ticks.length = unit(0.05, "cm"), # adjust length of tick (outside panel)
                axis.ticks =  element_line(linewidth = border.width)) # adjust width of tick(outside panel) (mm)
  
  return(p)
}

# RUN
p.scatter = lapply(pheno, FUN = perPhenoPlotScatterWithCommon,
                   avg.cpm.mtx = consensus.mtx.avg,
                   common.peaks,
                   common.bool = F,
                   scatter.fills,
                   logtf = F,
                   pt.size = 0.3,
                   keep.legend = F)



get_plot_limits <- function(plot) { # https://stackoverflow.com/a/40304848
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

axes.lims = lapply(p.scatter, FUN = get_plot_limits)  
ymax = lapply(axes.lims, FUN = function(x){x$ymax}) %>% unlist() %>% max()
xmax = lapply(axes.lims, FUN = function(x){x$xmax}) %>% unlist() %>% max()
ymin = lapply(axes.lims, FUN = function(x){x$ymin}) %>% unlist() %>% min()
xmin = lapply(axes.lims, FUN = function(x){x$xmin}) %>% unlist() %>% min()

# Save Plots ####

p.scatter.egg = lapply(p.scatter, FUN = set_panel_size, width = unit(1.5, units = "in"), 
                       height = unit(1.5, units = "in"))

for(i in 1:4){
  ggsave(filename = paste0(scratch.dir, "/figures/signalscatter_",pheno[i],".pdf"),
         plot = p.scatter[[i]] + 
           coord_cartesian(ylim = c(10^ymin,10^ymax),
                           xlim = c(10^xmin,10^xmax)),
         width = 2,
         height = 2,
         units = "in",
         device = "pdf")
}
