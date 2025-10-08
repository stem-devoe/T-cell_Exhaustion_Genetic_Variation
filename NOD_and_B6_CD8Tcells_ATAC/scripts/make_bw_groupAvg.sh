# CPM normalized and scaled to reciprocal FRIP

set -eu

samples=/Scottbrowne/members/smd/Projects/SD029/sample_lists/low-noise_archr-iterative-merge_input.tsv
blacklist=/Scottbrowne/members/smd/genomes/mm10/blacklist/mm10-blacklist-v2_chrUnRandomBlacklist.bed

bigwig_path=/scratch2/devoes/SD029/bigwig_cpm_recip-frip/ # make sure ends in "/"
bigwig_ext="_cpm_recip-frip.bw"

outdir=/scratch2/devoes/SD029/bigwig_cpm_recip-frip_groupAvg

mkdir -p $outdir


# get groups 
all_groups=$(cut -f2 $samples | sort | uniq | tr "\n" " ")
all_groups=${all_groups/" Group"/} # remove header

# for each group, make avg bw
for group in $all_groups
do

group_samples=$(grep $group $samples | cut -f1 | awk -v path="${bigwig_path}" -v ext="${bigwig_ext}" '{print path$1ext}' | tr "\n" " ") 

echo $group
echo $group_samples

bigwigAverage -b $group_samples \
-o ${outdir}/${group}_cpm_recip-frip_groupAvg.bw \
--binSize 10 \
--blackListFileName $blacklist \
--numberOfProcessors 8

done



