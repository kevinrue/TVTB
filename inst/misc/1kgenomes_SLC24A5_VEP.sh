#!/bin/sh

# In the /extdata folder of the TVTB package

# Download a 1000-genomes VCF file and index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

# Subset using genomic coordinates
bcftools view -r 15:48413169-48434869 -o SLC24A5.vcf ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Compress and index the VCF file
bgzip SLC24A5.vcf
tabix SLC24A5.vcf.gz

# Run the VEP script
variant_effect_predictor.pl \
    --cache \
    --offline \
    --everything \
    --fork 2 \
    --assembly GRCh37 \
    --plugin CADD \
    -i SLC24A5.vcf.gz \
    --vcf \
    -o chr15.phase3_integrated.vcf

# Compress and index the VCF file
bgzip chr15.phase3_integrated.vcf
tabix chr15.phase3_integrated.vcf.gz
