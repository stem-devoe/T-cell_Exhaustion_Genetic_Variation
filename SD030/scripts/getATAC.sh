# download atacseq fastq files
mkdir -p /Scottbrowne/seq/tmp/devoes/SD030/fastq/atac
cd /Scottbrowne/seq/tmp/devoes/SD030/fastq/atac

while read -r accession
do
echo $accession
fasterq-dump --split-files $accession -e 4 -O /Scottbrowne/seq/tmp/devoes/SD030/fastq/atac
done < /Scottbrowne/members/smd/Projects/SD030/SRA_Accession_ATAC.txt

gzip /Scottbrowne/seq/tmp/devoes/SD030/fastq/atac/*.fastq