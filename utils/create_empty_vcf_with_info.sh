#!/bin/bash
# create_empty_vcf_with_info.sh
# Usage: ./create_empty_vcf_with_info.sh snp_list.txt sample_list.txt output.vcf

snp_list="$1"      # file with SNPs: CHROM POS ID REF ALT
sample_list="$2"   # file with one sample name per line
out_vcf="$3"       # output VCF filename
bcftools="$4"

# 1. Create VCF header
vcf_header="##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">
##INFO=<ID=EAF,Number=A,Type=Float,Description=\"Effect allele frequency\">
##INFO=<ID=INFO_SCORE,Number=1,Type=Float,Description=\"Imputation info score\">
##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Hardy-Weinberg equilibrium p-value\">
##INFO=<ID=ERC,Number=A,Type=Integer,Description=\"Effect allele count\">
##INFO=<ID=EAC,Number=A,Type=Integer,Description=\"Effect allele count (alternate)\">
##INFO=<ID=PAF,Number=A,Type=Float,Description=\"Population allele frequency\">
##INFO=<ID=REF_PANEL,Number=1,Type=String,Description=\"Reference panel\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

# Append sample names
while read sample; do
    vcf_header="${vcf_header}\t${sample}"
done < "$sample_list"

# Write header to output VCF
echo -e "$vcf_header" > "$out_vcf"

# 2. Create body with empty genotypes and INFO fields
# SNP list format: CHROM POS ID REF ALT
while read chrom pos id ref alt; do
    info="EAF=.;INFO_SCORE=.;HWE=.;ERC=.;EAC=.;PAF=.;REF_PANEL=."
    line="${chrom}\t${pos}\t${id}\t${ref}\t${alt}\t.\t.\t${info}\tGT:GP:DS"
    while read sample; do
        line="${line}\t./.:.:."
    done < "$sample_list"
    echo -e "$line" >> "$out_vcf"
done < "$snp_list"

# 3. bgzip-compress and index with bcftools
bgzip -f "$out_vcf"
bcftools index -f "${out_vcf}.gz"

echo "Empty VCF with INFO fields created, gzipped, and indexed: ${out_vcf}.gz"