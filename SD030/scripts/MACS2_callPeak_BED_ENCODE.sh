set -eu

## Calling method implemented by the ENCODE pipeline
# https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_macs2_atac.py

# set paths 
chrSizes=/Scottbrowne/members/smd/genomes/mm10/mm10.chrom.sizes
bedPaths=/Scottbrowne/members/smd/Projects/SD030/sample_inputs/bed-paths_for-peak-call.txt
outdir=/Scottbrowne/seq/tmp/devoes/SD030/macs2_bed/peaks


while read -r bed
do

bedBase=$(basename $bed .tn5.tagAlign.gzip)
echo $bedBase

macs2 callpeak -t $bed -f BED --outdir $outdir -n $bedBase -g mm -q 0.01 --shift -75 --extsize 150 --nomodel --nolambda --keep-dup all --call-summits

done < $bedPaths


