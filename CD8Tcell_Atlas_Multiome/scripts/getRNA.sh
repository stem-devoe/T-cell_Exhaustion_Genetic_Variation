# download atacseq fastq files
mkdir -p /Scottbrowne/seq/tmp/devoes/SD030/fastq/rna
cd /Scottbrowne/seq/tmp/devoes/SD030/fastq/rna

while read -r accession
do
echo $accession
fasterq-dump --split-files $accession -e 4 -O /Scottbrowne/seq/tmp/devoes/SD030/fastq/rna
done < /Scottbrowne/members/smd/Projects/SD030/SRA_Accession_RNA.txt

gzip /Scottbrowne/seq/tmp/devoes/SD030/fastq/rna/*.fastq