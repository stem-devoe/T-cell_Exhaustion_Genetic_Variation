set -eux
# GV sig peak intersection with GRANIE p2g

save_date="2024-10-06"

# granie peak to gene correlations (not filtered by rho)
p2g=/Scottbrowne/members/smd/Projects/SD030/data/GRANIE/${save_date}_peak2gene_pearson.bed
# path to GV sig peaks
sigPath=/Scottbrowne/members/smd/Projects/SD029/data/150bp/padj_0.1/sig_peaks/FC_1.5
#
acqPath=/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/150bp/padj_0.1/FC_1.5
#SNPs bed files
nodAllSNP=/Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/vcf/bed/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.PASS.hqHetSNP.vcf.bed
nodHetSNP=/Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/vcf/bed/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.hqHetSNP.vcf.bed
nodOnlyPassSNP=/Scottbrowne/members/smd/genomes/NOD_SHILTJ_v5_REL1505/vcf/bed/NOD_ShiLtJ.mgp.v5.snps.dbSNP142.sort.chrUCSC.onlyPASS.vcf.bed
# high density HetSNP regions
#hdHetSNP=/scratch2/devoes/SD029/hSNP/dense-hSNP_from_mm10_5kb-window_1kb-step.merged.bed
hdHetSNP=/Scottbrowne/members/smd/Projects/SD032/dense-hSNP_from_mm10_5kb-window_1kb-step.merged.bed
# strain specific motif sites
B6SpecificTFBS=/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_consensus_b6sites.bed
NODSpecificTFBS=/Scottbrowne/members/smd/Projects/SD029/motif_disruption/HOCOMOCO/hocomocov11_consensus_nodsites.bed

outPath=/Scottbrowne/members/smd/Projects/SD029/GRANIE_p2g_intersections

mkdir -p $outPath


for sigBed in ${sigPath}/NOD-B6*.up.bed
do
sampleOut=${outPath}/$(basename $sigBed .bed)
bedtools intersect -wo -a $sigBed -b $p2g > ${sampleOut}.GRANIE-p2g.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodAllSNP > ${sampleOut}.GRANIE-p2g.AllSNP.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodHetSNP > ${sampleOut}.GRANIE-p2g.HetSNP.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodOnlyPassSNP > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $B6SpecificTFBS > ${sampleOut}.GRANIE-p2g.B6SpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $NODSpecificTFBS > ${sampleOut}.GRANIE-p2g.NODSpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed -b $B6SpecificTFBS > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.B6SpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed -b $NODSpecificTFBS > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.NODSpecificTFBS.bed
# below does not guarantee the peak to have a polymorphism
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $hdHetSNP > ${sampleOut}.GRANIE-p2g.hdHetSNP.bed
done

for sigBed in ${sigPath}/NOD-B6*.down.bed
do
sampleOut=${outPath}/$(basename $sigBed .bed)
bedtools intersect -wo -a $sigBed -b $p2g > ${sampleOut}.GRANIE-p2g.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodAllSNP > ${sampleOut}.GRANIE-p2g.AllSNP.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodHetSNP > ${sampleOut}.GRANIE-p2g.HetSNP.bed
bedtools intersect -wa -u -a ${sampleOut}.GRANIE-p2g.bed -b $nodOnlyPassSNP > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $B6SpecificTFBS > ${sampleOut}.GRANIE-p2g.B6SpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $NODSpecificTFBS > ${sampleOut}.GRANIE-p2g.NODSpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed -b $B6SpecificTFBS > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.B6SpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.OnlyPassSNP.bed -b $NODSpecificTFBS > ${sampleOut}.GRANIE-p2g.OnlyPassSNP.NODSpecificTFBS.bed
# below does not guarantee the peak to have a polymorphism
bedtools intersect -wo -a ${sampleOut}.GRANIE-p2g.bed -b $hdHetSNP > ${sampleOut}.GRANIE-p2g.hdHetSNP.bed
done

# for peaks acquired in exhausted state only
for i in LCMVcl13 CART
do
acqBed=${acqPath}/${i/"cl13"/}/Gained_${i/"cl13"/}-to-Naive_any-strain_peak_subset.bed # boo on me for inconsistent nomenclature
sampleOut=${outPath}/acquired
mkdir -p $sampleOut
for j in up down
do
bedtools intersect -wo -a $acqBed -b $p2g | bedtools intersect -wa -a - -b ${sigPath}/NOD-B6_${i}.${j}.bed > ${sampleOut}/acquired-${i}.GRANIE-p2g.NOD-B6-${j}.bed
bedtools intersect -wo -a ${sampleOut}/acquired-${i}.GRANIE-p2g.NOD-B6-${j}.bed -b $B6SpecificTFBS | bedtools intersect -wa -a - -b $nodOnlyPassSNP > ${sampleOut}/acquired-${i}.GRANIE-p2g.NOD-B6-${j}.OnlyPassSNP.B6SpecificTFBS.bed
bedtools intersect -wo -a ${sampleOut}/acquired-${i}.GRANIE-p2g.NOD-B6-${j}.bed -b $NODSpecificTFBS | bedtools intersect -wa -a - -b $nodOnlyPassSNP > ${sampleOut}/acquired-${i}.GRANIE-p2g.NOD-B6-${j}.OnlyPassSNP.NODSpecificTFBS.bed
done
done
