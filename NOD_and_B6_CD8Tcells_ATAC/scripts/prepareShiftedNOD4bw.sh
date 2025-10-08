for i in /scratch2/devoes/SD029/nfcore_outs/NOD/shifted/bowtie2/*SHILTJ.bam
do
	samtools sort $i -o ${i/"sorted_shifted_from_NOD_SHILTJ.bam"/"shifted_from_NOD_SHILTJ.coordinate-sorted.bam"}
	samtools index ${i/"sorted_shifted_from_NOD_SHILTJ.bam"/"shifted_from_NOD_SHILTJ.coordinate-sorted.bam"}
done
