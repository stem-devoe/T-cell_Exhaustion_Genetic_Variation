set -eu

## run from outdir because MMARGE will place files in the PWD I think

# define paths and files
samples=/Scottbrowne/members/smd/Projects/SD029/sample_lists/NOD_samples.txt
outdir=/scratch2/devoes/SD029/nfcore_outs/NOD/shifted/bowtie2
nfcoreDir=/scratch2/devoes/SD029/nfcore_outs/NOD/bowtie2/merged_library

# enable activation of conda environments
source /home/devoes/miniconda3/etc/profile.d/conda.sh
conda activate CUTnTAG

# make output directory 
mkdir -p $outdir

while read -r lines
do

# strip fasta path of sample
lines=$(basename $lines)
echo $lines

# set bam name
bam=${nfcoreDir}/${lines}_REP1.mLb.clN.sorted.bam
echo $bam

# save header and remove "_allele_1"
# shifting with header strips after all "_" and messes up the nonstandard chr in the header when we try to convert back to bam
samtools view -H $bam > ${outdir}/${lines}.samheader
sed -i 's!_allele_1!!g' ${outdir}/${lines}.samheader

# convert to sam
sam=${outdir}/${lines}_REP1.mLb.clN.sorted.sam
samtools view $bam -o $sam
echo $sam

#shift using MMARGE back to B6 coordinates
/Scottbrowne/members/smd/MMARGE_v1.0/bin/MMARGE.pl shift -sam -ind NOD_SHILTJ -data_dir /Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/MMARGE/mutation_files -files $sam

# add header back onto shifted sam file and convert sam to bam
cat ${outdir}/${lines}.samheader ${sam/".sam"/"_shifted_from_NOD_SHILTJ.sam"} | samtools view -hb - -o ${sam/".sam"/"_shifted_from_NOD_SHILTJ.bam"}

# remove sam files and sam header file
rm $sam
rm ${sam/".sam"/"_shifted_from_NOD_SHILTJ.sam"}
rm ${outdir}/${lines}.samheader

done < $samples



