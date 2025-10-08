set -eu
## Calling method implemented by the ENCODE pipeline
# https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_bam2ta.py#L95


# set paths 
bamPaths=/Scottbrowne/members/smd/Projects/SD030/sample_inputs/bam-paths_tn5-shift_2bed.txt
outdir=/Scottbrowne/seq/tmp/devoes/SD030/macs2_bed/bed

mkdir -p $outdir

while read -r bam
do

bamBase=$(basename $bam .bam) # B6 files are .sorted.bam
bamBase=${bamBase/'.sorted'/}
echo $bamBase

## get name sorted bam
samtools sort -O bam -n $bam -o ${outdir}/${bamBase}.namesort.bam 

## convert bampe to bedpe
bedtools bamtobed -bedpe -mate1 -i ${outdir}/${bamBase}.namesort.bam  | gzip -nc > ${outdir}/${bamBase}.bedpe.gzip

rm ${outdir}/${bamBase}.namesort.bam

## give each read its own bed entry
# places "N" as placeholder in name field
# places 1000 into score field
# re-did the ENCODE syntax so it's easier to read: awk 'BEGIN{OFS="\t"}{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}'

zcat -f ${outdir}/${bamBase}.bedpe.gzip | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N",1000,$9"\n"$4,$5,$6,"N",1000,$10}' | \
gzip -nc > ${outdir}/${bamBase}.tagAlign.gzip


## shift for tn5 insertion and format bed
# postitive strand reads get +4
# negative strand reads get -5 
# check that reads are in proper orientation (ie start is lowest genomic coordinate and end is further downstream)
#     if not, set reads to 1 bp
#         set pos strand start to end - 1 (because end is inclusive but start is exclusive in bam. bam is 0-based and half-open)
#         and set neg strand end to start + 1 (to make the new end coordinate inclusive via bed0-based half-open format)

zcat -f ${outdir}/${bamBase}.tagAlign.gzip | \
awk 'BEGIN{OFS="\t"}{if($6 == "+"){$2 = $2 +4} else if ($6 == "-"){$3 = $3 - 5}; if($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}} print}' | \
gzip -nc > ${outdir}/${bamBase}.tn5.tagAlign.gzip


done < $bamPaths
