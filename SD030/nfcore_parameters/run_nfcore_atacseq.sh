nextflow run nf-core/atacseq -profile docker \
-r 2.1.2 \
-c /mnt/nautilus/members/smd/Projects/SD030/nfcore_parameters/nfcore_atacseq_custom.config \
-params-file /mnt/nautilus/members/smd/Projects/SD030/nfcore_parameters/nfcore_atacseq.yaml \
--input /mnt/nautilus/members/smd/Projects/SD030/nfcore_parameters/atac_nfcore_input.csv \
--outdir /mnt/nautilus/seq/tmp/devoes/SD030/nfcore_outs/atac