
set -eux

exonIDFile=/Scottbrowne/seq/tmp/devoes/SD032/Ccl27_Ccl21_Il11ra/exonIDs.txt
alignmentPerc=0.8
blastDB=GCA_921998325.2_NOD_ShiLtJ_v3_genomic
blastDBPath=/Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJEB47108/blast_db/GCA_921998325.2_NOD_ShiLtJ_v3_genomic
exonFasta=/Scottbrowne/seq/tmp/devoes/SD032/Ccl27_Ccl21_Il11ra/exons.fa
chrSizes=/Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJEB47108/GCA_921998325.2_NOD_ShiLtJ_v3/GCA_921998325.2_NOD_ShiLtJ_v3_genomic.sizes
outPathBase=/Scottbrowne/seq/tmp/devoes/SD032/Ccl27_Ccl21_Il11ra/exons_blast/NOD

mkdir -p $outPathBase

outCopies=${outPathBase}/copy_coords.tsv
copyCounts=${outPathBase}/copy_counts.tsv

: > $outCopies


while read -r exonID 
do

outPath=${outPathBase}/${exonID}
mkdir -p $outPath

cat $exonFasta | grep $exonID -A 1 | blastn -task megablast -db $blastDBPath -query - -out ${outPath}/${exonID}-megablastHits-${blastDB}.tsv -evalue 0.05  -outfmt '6 saccver sstart send evalue score length sstrand qlen pident nident mismatch gapopen gaps sseq qframe sframe' -num_threads 4 

# output to sorted bed and filter hits by alignment length
cut -f 1,2,3,6,7,8 ${outPath}/${exonID}-megablastHits-${blastDB}.tsv | sort -k1,1 -k 2,2n | awk -v alignPerc="$alignmentPerc" -v exonID="$exonID" 'OFS=FS="\t" {if(($4/$6) >= alignPerc) print $1, $2 < $3 ? $2 : $3, $2 < $3 ? $3 : $2,exonID,$5}' > ${outPath}/${exonID}-megablastHits-${blastDB}.bed

awk 'OFS=FS="\t" {print $1,$2,$3,$4,$5,FNR}' ${outPath}/${exonID}-megablastHits-${blastDB}.bed >> $outCopies

done < $exonIDFile

bedtools sort -g $chrSizes -i $outCopies > ${outCopies/".tsv"/".bed"}
