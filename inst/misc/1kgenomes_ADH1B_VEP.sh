#!/bin/sh

# In the /extdata folder of the TVTB package

# Download a 1000-genomes VCF file and index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

# Subset using genomic coordinates
bcftools view -r 4:100226121-100242558 -o ADH1B.vcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Compress and index the VCF file
bgzip ADH1B.vcf
tabix ADH1B.vcf.gz

# Run the VEP script
variant_effect_predictor.pl \
    --cache \
    --offline \
    --everything \
    --fork 2 \
    --assembly GRCh37 \
    --plugin CADD \
    -i ADH1B.vcf.gz \
    --vcf \
    -o chr4.phase3_integrated.vcf

# Compress and index the VCF file
bgzip chr4.phase3_integrated.vcf
tabix chr4.phase3_integrated.vcf.gz

# When the process is complete:
# Manually remove ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz*
