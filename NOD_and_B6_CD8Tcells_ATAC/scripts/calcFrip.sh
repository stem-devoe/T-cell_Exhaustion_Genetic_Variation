set -eu

## determine FRiP scores for each sample

# nfcore bamfile sample name through '_REP1.mLb.clN' or '_REP1.mLb.clN.shifted_from_NOD_SHILTJ'
sampleNames=/Scottbrowne/members/smd/Projects/SD029/sample_lists/samples_frip-scoring.txt
readPath=/scratch2/devoes/SD029/macs2_bed/bed # do we want reads or tn5 shifted reads to count (tagAlign or tn5_tagAlign folder)
peaks=/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.bed
outFile=/Scottbrowne/members/smd/Projects/SD029/data/FRIP_ConsensusPeaks-150bpExt.txt


while read -r sample
do

inters=$(bedtools intersect -a $readPath/${sample}.tagAlign.gzip -b $peaks -wa -u | wc -l)
reads=$(zcat $readPath/${sample}.tagAlign.gzip | wc -l)
frip=$(awk -v var1="${inters}" -v var2="${reads}" 'BEGIN {print var1/var2}')  # since there is no input file, need to tell awk where to begin

echo ${sample}$':\t'${frip} >> $outFile

done < $sampleNames
