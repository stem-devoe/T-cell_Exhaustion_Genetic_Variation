# should I add subject strand to output? Would this help with IDing inversions?

set -eux

peakList=/Scottbrowne/seq/tmp/devoes/SD032/Cd244_CopyBoundaries/peakIDs.txt
ArchRExt=150
blastDB=GCA_921998325.2_NOD_ShiLtJ_v3_genomic
blastDBPath=/Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJEB47108/blast_db/GCA_921998325.2_NOD_ShiLtJ_v3_genomic
peakFasta=/scratch2/devoes/SD029/macs2_bed/archr-iterative-merge_peaks-by-group_variableExt/globallambda/150bp/All_Samples.fwp.filter.non_overlapping.renamed.fa
outPathBase=/Scottbrowne/seq/tmp/devoes/SD032/Cd244_CopyBoundaries_Strand

outCopies=${outPathBase}/copy_coords.tsv

: > $outCopies

while read -r peakID
do

	outPath=${outPathBase}/${peakID}
	mkdir -p $outPath

	cat $peakFasta | grep $peakID -A 1 | blastn -task megablast -db $blastDBPath -query - -out ${outPath}/${peakID}-megablastHits-${blastDB}.tsv -evalue 0.05  -outfmt '6 saccver sstart send evalue score length sstrand pident nident mismatch gapopen gaps sseq qframe sframe' -num_threads 4 

	# output to sorted bed and filter hits by alignment length
	cut -f 1,2,3,6,7 ${outPath}/${peakID}-megablastHits-${blastDB}.tsv | sort -k1,1 -k 2,2n | awk -v ArchRExt="$ArchRExt" -v peakID="$peakID" 'OFS=FS="\t" {if($4 >= ((2 * ArchRExt) + 1) * 0.8) print $1, $2 < $3 ? $2 : $3, $2 < $3 ? $3 : $2,peakID,$5}' > ${outPath}/${peakID}-megablastHits-${blastDB}.bed

	awk 'OFS=FS="\t" {print $1,$2,$3,$4,$5,FNR}' ${outPath}/${peakID}-megablastHits-${blastDB}.bed >> $outCopies

done < $peakList
