nextflow run nf-core/atacseq -profile docker \
        -r 2.1.2 \
        -c /Scottbrowne/members/smd/Projects/SD032/nfcore_parameters/atacseq_custom.config \
        -params-file /Scottbrowne/members/smd/Projects/SD032/nfcore_parameters/nfcore_atacseq.yaml \
        --input /Scottbrowne/members/smd/Projects/SD032/nfcore_parameters/nfcore_atacseq-input.csv \
        --outdir /Scottbrowne/seq/tmp/devoes/SD032/nfcore_outs/GCA_921998325.2_NOD_ShiLtJ_v3_genomic \
		--max_cpus 24 \
		--max_memory 160.GB