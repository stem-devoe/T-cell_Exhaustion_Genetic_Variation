nextflow run nf-core/atacseq -profile docker \
        -r 2.1.2 \
        -c /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/atacseq_custom.config \
        -params-file /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_NOD.yaml \
        --input /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_NOD-input.csv \
        --outdir /scratch2/devoes/SD029/nfcore_outs/NOD \
	--max_cpus 30 \
	--max_memory 200.GB

nextflow run nf-core/atacseq -profile docker \
        -r 2.1.2 \
        -c /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/atacseq_custom.config \
        -params-file /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_B6.yaml \
        --input /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_B6-input.csv \
        --outdir /scratch2/devoes/SD029/nfcore_outs/B6 \
	--max_cpus 30 \
	--max_memory 200.GB

nextflow run nf-core/atacseq -profile docker \
        -r 2.1.2 \
        -c /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/atacseq_custom.config \
        -params-file /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_B6.yaml \
        --input /Scottbrowne/members/smd/Projects/SD029/nfcore_parameters/nfcore_atacseq_NOD-input.csv \
        --outdir /scratch2/devoes/SD029/nfcore_outs/NOD-to-B6 \
	--max_cpus 30 \
	--max_memory 200.GB
		

