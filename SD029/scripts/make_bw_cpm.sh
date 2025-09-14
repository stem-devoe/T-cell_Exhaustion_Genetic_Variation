set -eu

# shifted NOD bam files need to be coordinate sorted and indexed before creating bigwigs
# see prepareShiftedNOD4bw.sh

bamPathB6=/scratch2/devoes/SD029/nfcore_outs/B6/bowtie2/merged_library
bamPathNOD=/scratch2/devoes/SD029/nfcore_outs/NOD/shifted/bowtie2
outPath=/scratch2/devoes/SD029
sampleFRIP=/Scottbrowne/members/smd/Projects/SD029/sample_lists/bw-frip_input.csv # nfcore bamfile base through REP1.mlB.clN(.shifted_from_NOD_SHILTJ)
blacklist=/Scottbrowne/members/smd/genomes/mm10/blacklist/mm10-blacklist-v2_chrUnRandomBlacklist.bed

mkdir -p ${outPath}/bigwig_cpm
mkdir -p ${outPath}/bigwig_cpm_recip-frip

while read -r lines
do

# input is format of sampleName,FRIP,RecriprocalFRIP
sample=${lines%%,*}
recfrip=${lines##*,}

if [[ $sample == *"shifted"* ]]
then
bam=${bamPathNOD}/${sample}.coordinate-sorted.bam
else
bam=${bamPathB6}/${sample}.sorted.bam
fi

out_cpm=${outPath}/bigwig_cpm/${sample}_cpm.bw
out_cpm_rf=${outPath}/bigwig_cpm_recip-frip/${sample}_cpm_recip-frip.bw

echo $bam
echo $recfrip
echo $out_cpm
echo $out_cpm_rf

bamCoverage -b $bam \
-o $out_cpm \
--binSize 10 \
--blackListFileName $blacklist \
--numberOfProcessors 8 \
--normalizeUsing CPM

bamCoverage -b $bam \
-o $out_cpm_rf \
--binSize 10 \
--blackListFileName $blacklist \
--numberOfProcessors 8 \
--normalizeUsing CPM \
--scaleFactor $recfrip

done < $sampleFRIP