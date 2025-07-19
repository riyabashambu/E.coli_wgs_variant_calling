# E.coli_wgs_variant_calling
The purpose of this project is to identify genetic variations, specifically Single Nucleotide Polymorphisms (SNPs) and insertions/deletions (indels), in Escherichia coli strains using whole genome sequencing (WGS) data. This analysis helps in understanding strain-specific mutations, evolutionary differences, and functional changes that may influence pathogenicity, drug resistance, or adaptation.

By implementing a standardized bioinformatics pipeline using tools like BWA, SAMtools, GATK, and IGV, this project demonstrates the end-to-end process of:

-Quality control of raw sequence reads

-Reference-based genome alignment

-Variant calling and filtering

-Visualization and interpretation of variants

This project can serve as a template or learning resource for bacterial genome variant analysis using short-read sequencing data.

*RECOMMENDED ENVIRONMENT SETUP*

 we can install all required tools using Conda:
 
'''conda create -n ecoli_variant_env \
  fastqc trimmomatic bwa samtools gatk4 snpeff bcftools igv \
  -c bioconda -c conda-forge -y

conda activate ecoli_variant_env'''


