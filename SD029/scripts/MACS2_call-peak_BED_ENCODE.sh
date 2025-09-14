set -eu

## Calling method implemented by the ENCODE pipeline
# https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_macs2_atac.py

# set paths 
chrSizes=/Scottbrowne/members/smd/genomes/mm10/mm10.chrom.sizes
bedPaths=/Scottbrowne/members/smd/Projects/SD029/sample_lists/bed-paths_for-peak-call.txt
outdir=/scratch2/devoes/SD029/macs2_bed/peaks

mkdir -p $outdir

while read -r bed
do

bedBase=$(basename $bed .tn5.tagAlign.gzip)
echo $bedBase

macs2 callpeak -t $bed -f BED --outdir $outdir -n $bedBase -g mm -q 0.01 --shift -75 --extsize 150 --nomodel --nolambda --keep-dup all --call-summits

done < $bedPaths


