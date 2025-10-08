nextflow run nf-core/rnaseq -profile docker \
-r 3.13.2 \
-params-file /mnt/nautilus/members/smd/Projects/SD030/nfcore_parameters/nfcore_rnaseq.yaml \
--input /mnt/nautilus/members/smd/Projects/SD030/nfcore_parameters/rna_nfcore_input.csv \
--outdir /mnt/nautilus/seq/tmp/devoes/SD030/nfcore_outs/rna