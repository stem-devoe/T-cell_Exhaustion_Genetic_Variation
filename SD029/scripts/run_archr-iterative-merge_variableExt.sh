# https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging

# "score-per-million cutoff to be used in determining which peak calls represent true peaks. The higher the spm cutoff, the more significant a peak call has to be to be retained. We recommended setting spm to between 2 and 5 as a starting point. Default value is 5"



spm_min=4 # i used a q value cutoff of 0.01 instead of p value (as recommended in ArchR docs) so I am going to lower the stringency for spm (default 5) - the score value is based on int(-10*log10pvalue) or int(-10*log10qvalue) depending on the callpeak filtering stat you used: p or q

for extSize in 50 100 125 150 200 250 500
do
	mkdir -p /scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/${extSize}bp

	Rscript /Scottbrowne/members/smd/Projects/SD022/scripts/createIterativeOverlapPeakSet.R \
		--metadata /Scottbrowne/members/smd/Projects/SD029/sample_lists/low-noise_archr-iterative-merge_input.tsv \
		--macs2dir /scratch2/devoes/SD029/macs2_bed/peaks/ \
		--outdir /scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/${extSize}bp/ \
		--suffix _summits.bed \
		--blacklist /Scottbrowne/members/smd/genomes/mm10/blacklist/mm10-blacklist-v2_chrUnRandomBlacklist.bed \
		--genome mm10 \
		--spm $spm_min \
		--rule "(n+1)/2" \
		--extend $extSize

done
