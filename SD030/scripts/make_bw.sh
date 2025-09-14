# cpm normalize and reciprocal FRIP scale

set -eu
binSize=10
bamPath=/Scottbrowne/seq/tmp/SD030/nfcore_outs/atac/bowtie2/merged_library
outPath=/Scottbrowne/seq/tmp/devoes/SD030/bigwig_binSize${binSize}/bigwig_cpm_recip-frip

# format for input is: sampleName,FRIP,ReciprocalFRIP
sampleFRIP=/Scottbrowne/members/smd/Projects/SD030/sample_inputs/bw-frip_input.csv # sampleName is nfcore bamfile base through REP1.mlB.clN

blacklist=/Scottbrowne/members/smd/genomes/mm10/blacklist/mm10-blacklist-v2_chrUnRandomBlacklist.bed

mkdir -p $outPath

while read -r lines
do

# input is format of sampleName,FRIP,RecriprocalFRIP
sample=${lines%%,*}
recfrip=${lines##*,} 

bam=${bamPath}/${sample}.sorted.bam

out=${outPath}/${sample}_cpm_recip-frip.bw

echo $bam
echo $recfrip
echo $out

bamCoverage -b $bam \
-o $out \
--binSize ${binSize} \
--blackListFileName $blacklist \
--numberOfProcessors 8 \
--normalizeUsing CPM \
--scaleFactor $recfrip

done < $sampleFRIP

