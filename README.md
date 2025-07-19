# E.coli_wgs_variant_calling
The purpose of this project is to identify genetic variations, specifically Single Nucleotide Polymorphisms (SNPs) and insertions/deletions (indels), in Escherichia coli strains using whole genome sequencing (WGS) data. This analysis helps in understanding strain-specific mutations, evolutionary differences, and functional changes that may influence pathogenicity, drug resistance, or adaptation.

By implementing a standardized bioinformatics pipeline using tools like BWA, SAMtools, GATK, and IGV, this project demonstrates the end-to-end process of:
1. Quality control of raw sequence reads 
2. Reference-based genome alignment 
3. Variant calling and filtering 
4. Visualization and interpretation of variants 

This project can serve as a template or learning resource for bacterial genome variant analysis using short-read sequencing data.


**RECOMMENDED ENVIRONMENT SETUP**

 we can install all required tools using Conda:
 
<pre>bash
 
 conda create -n ecoli_variant_env \
  fastqc trimmomatic bwa samtools gatk4 snpeff bcftools igv \
  -c bioconda -c conda-forge -y

conda activate ecoli_variant_env</pre>


**TOOLS USED**

1. FastQC 
2. Trimmomatic 
3. BWA 
4. SAMtools 
5. GATK 
6. SnpEff or BCFtools annotate 
7. IGV 


**FILE FORMT USED**

| **Format** | **Extension**       | **Description**                                           |
|------------|---------------------|-----------------------------------------------------------|
| FASTQ      | `.fastq`, `.fastq.gz` | Raw sequencing reads (e.g., paired-end: R1 and R2)      |
| FASTA      | `.fasta`, `.fa`     | Reference genome sequence                                 |
| SAM        | `.sam`              | Alignment file from BWA (text format)                     |
| BAM        | `.bam`              | Binary alignment file (compressed SAM)                    |
| BAI        | `.bai`              | Index file for BAM (used for visualization in IGV)        |
| VCF        | `.vcf`              | Variant Call Format – lists SNPs and indels               |
| VCF.GZ     | `.vcf.gz`           | Compressed VCF file                                       |
| GFF / GTF  | `.gff`, `.gtf`      | Genome annotation file (used for variant annotation)      |
| HTML       | `.html`             | Quality control reports from FastQC                       |
| TXT        | `.txt`              | Logs, notes, or summary output (optional)                 |


**WORKFLOW**

The overview of the workflow is-
1. Quality check of raw reads using FastQC  
2. Read trimming with Trimmomatic  
3. Alignment to reference genome using BWA  
4. SAM to BAM conversion and sorting using SAMtools  
5. Variant calling using GATK  
6. Variant annotation using SnpEff  
7. Visualization in IGV  

**RAW DATA AND REFFERENCE GENOME**

Raw Sequencing Data

 Organism: *Escherichia coli* (strain K-12 substr. MG1655) 
 Source: European Nucleotide Archive (ENA) 
 Run Accession: [SRR3191544](https://www.ebi.ac.uk/ena/browser/view/SRR3191544) 
 Sequencing Platform: Illumina HiSeq  
 Read Type: Paired-end (2 × 150 bp) 

<pre>bash
# Read 1
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/004/SRR3191544/SRR3191544_1.fastq.gz

# Read 2
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/004/SRR3191544/SRR3191544_2.fastq.gz </pre>

Reference Genome

  Reference Strain: E. coli K-12 substrain MG1655 
  Source: NCBI RefSeq
  Assembly Accession: NC_000913.3
  File: ecoli_reference.fasta

<pre>bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_reference.fasta </pre>
