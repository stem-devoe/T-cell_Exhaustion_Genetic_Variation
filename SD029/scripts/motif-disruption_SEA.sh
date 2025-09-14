# make sure meme/bin is in path: export PATH=/Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/SEA/meme-5.5.7/meme/bin:$PATH  
set -eux
for i in CART LCMV Naive
do
for j in ETS Runt bZIP # SP1-KLF
do
mkdir -p /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/SEA/${i}/${j} 
sea  --p /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/homer/${i}/${j}/target.fa --m /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/SEA/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme --n /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/homer/${i}/${j}/background.fa --notrim --oc /Scottbrowne/seq/tmp/devoes/SD029/motif_disruption/SEA/${i}/${j} --qvalue --thresh 1.0
done
done