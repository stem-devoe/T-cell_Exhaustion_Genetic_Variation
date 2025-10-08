## Build blast databases ##

zcat /Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJNA944543/GCA_031763685.1_NOD_SCID/GCA_031763685.1_NOD_SCID_genomic.fna.gz | makeblastdb -in - -parse_seqids -dbtype nucl -title GCA_031763685.1_NOD_SCID_genomic -blastdb_version 5 -logfile /Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJNA944543/blast_db/GCA_031763685.1_NOD_SCID_genomic.makeblastdb.log -out GCA_031763685.1_NOD_SCID_genomic

zcat /Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJEB47108/GCA_921998325.2_NOD_ShiLtJ_v3/GCA_921998325.2_NOD_ShiLtJ_v3_genomic.fna.gz | makeblastdb -in - -parse_seqids -dbtype nucl -title GCA_921998325.2_NOD_ShiLtJ_v3_genomic -blastdb_version 5 -logfile /Scottbrowne/seq/tmp/devoes/SD032/BioProject_PRJEB47108/blast_db/GCA_921998325.2_NOD_ShiLtJ_v3_genomic.makeblastdb.log -out GCA_921998325.2_NOD_ShiLtJ_v3_genomic

makeblastdb -in /Scottbrowne/members/smd/genomes/mm10/fasta/mm10.fa -parse_seqids -dbtype nucl -title mm10 -blastdb_version 5 -logfile /Scottbrowne/seq/tmp/devoes/SD032/mm10/blast_db/mm10.makeblastdb.log -out mm10


cat /Scottbrowne/seq/tmp/devoes/mm10.p6/chroms/*.fa | makeblastdb -in - -parse_seqids -dbtype nucl -title mm10.p6 -blastdb_version 5 -logfile /Scottbrowne/seq/tmp/devoes/mm10.p6/blast_db/mm10.p6.makeblastdb.log -out mm10.p6




