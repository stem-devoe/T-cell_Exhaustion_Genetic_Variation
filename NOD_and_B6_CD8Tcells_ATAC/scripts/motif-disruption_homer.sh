# conda activate homer
set -eux
date="2024-10-05"
for i in CART LCMV Naive
do
for j in ETS Runt bZIP # SP1-KLF
do
/home/devoes/miniconda3/envs/homer/bin/findMotifsGenome.pl /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/homer/${i}/${date}_${i}_${j}_peaks4homer.bed mm10 /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/homer/${i}/${j} -size 300 -mask -p 8 -dumpFasta
done
done