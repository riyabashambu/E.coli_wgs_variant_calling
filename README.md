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

**STEP 1. RAW DATA AND REFFERENCE GENOME**

Raw Sequencing Data

-Organism: *Escherichia coli* (strain K-12 substr. MG1655)  
-Source: European Nucleotide Archive (ENA)  
-Run Accession: [SRR3191544](https://www.ebi.ac.uk/ena/browser/view/SRR3191544)    
-Sequencing Platform: Illumina HiSeq   
 -Read Type: Paired-end (2 × 150 bp) 

<pre>bash
# Read 1
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/004/SRR3191544/SRR3191544_1.fastq.gz

# Read 2
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/004/SRR3191544/SRR3191544_2.fastq.gz </pre>

Reference Genome

 -Reference Strain: E. coli K-12 substrain MG1655 
 -Source: NCBI RefSeq
 -Assembly Accession: NC_000913.3
 -File: ecoli_reference.fasta

<pre>bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_reference.fasta </pre>

**STEP 2. QUALITY CONTROL OF RAW READS**
Objective:
To check the quality of sequencing reads to identify low-quality bases, adapter contamination, and any other issues before alignment.

Tool: FastQC

<pre>bash

# Run FastQC on both paired-end FASTQ files
fastqc SRR3191544_1.fastq.gz SRR3191544_2.fastq.gz </pre>

Output:
FastQC generates .html and .zip reports for each input file. Open the .html files in a browser to inspect:

-Per base quality scores 
-Adapter content
-GC content 
-Overrepresented sequences

Place output files into a dedicated directory:

<pre>bash

mkdir -p qc_reports
fastqc SRR3191544_1.fastq.gz SRR3191544_2.fastq.gz -o qc_reports/ </pre>

**STEP 3. TRIMMING LOW QUALITU READS AND ADAPTERS**

Objective: Remove sequencing adapters and trim low-quality bases from raw reads to improve downstream alignment accuracy.

Tool: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

 Input:
- `SRR3191544_1.fastq.gz` (Read 1)
- `SRR3191544_2.fastq.gz` (Read 2)

<pre> bash
trimmomatic PE -phred33 \
  SRR3191544_1.fastq.gz SRR3191544_2.fastq.gz \
  trimmed_R1_paired.fq.gz trimmed_R1_unpaired.fq.gz \
  trimmed_R2_paired.fq.gz trimmed_R2_unpaired.fq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  SLIDINGWINDOW:4:20 MINLEN:50 </pre>

Output:
-trimmed_R1_paired.fq.gz 
-trimmed_R2_paired.fq.gz 

**STEP 4. ALIGNING READS TO THE REFERENCE GENOME**

Objective: Map the sequencing reads to the *Escherichia coli* reference genome to determine their genomic positions.

Tool: [BWA-MEM](http://bio-bwa.sourceforge.net/) — a fast and accurate aligner for short reads.

Input:
- Reference genome: `ecoli_reference.fasta`
- Trimmed paired-end reads:
  - `trimmed_R1_paired.fq.gz`
  - `trimmed_R2_paired.fq.gz`

##1. Index the reference genome
<pre> bash
bwa index ecoli_reference.fasta </pre>
##2. Align paired-end reads
<pre>bash
bwa mem ecoli_reference.fasta \
  trimmed_R1_paired.fq.gz trimmed_R2_paired.fq.gz > aligned.sam </pre>

**STEP 5. SAM TO BAM CONVERSION, SORTING AND INDEXING**

Objective: Prepare the alignment file for variant calling by converting it to a compressed binary format (BAM), sorting the reads, and creating an index.

Tools: [SAMtools](http://www.htslib.org/)

Input: - `aligned.sam` — SAM file generated by BWA

##1. Convert SAM to BAM
<pre>bash
samtools view -S -b aligned.sam > aligned.bam </pre>

##2. Sort BAM file
<pre>bash
samtools sort aligned.bam -o aligned_sorted.bam </pre>

##3. Index the sorted BAM file
<pre>bash
samtools index aligned_sorted.bam </pre>

Output:
aligned_sorted.bam — sorted BAM file

aligned_sorted.bam.bai — BAM index file

**STEP 6. Variant Calling**

Objective: Identify SNPs and indels in the aligned *E. coli* genome using GATK's HaplotypeCaller.
Tool: [GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/)

🔧 Before running GATK, ensure to have:
> - A **sorted BAM** file with its **index**
> - A **reference genome** indexed for GATK:
>   - `.fasta`, `.fasta.fai`, `.dict`

##1. Create Sequence Dictionary for reference
<pre>bash
gatk CreateSequenceDictionary -R ecoli_reference.fasta </pre>

##2. Index the reference
<pre>bash
samtools faidx ecoli_reference.fasta </pre>

##3. Run HaplotypeCaller
<pre>bash
gatk HaplotypeCaller \
   -R ecoli_reference.fasta \
   -I aligned_sorted.bam \
   -O raw_variants.vcf \
   -ploidy 1 </pre>
Since E. coli is haploid, set -ploidy 1.

Output:
raw_variants.vcf — list of called variants (SNPs and indels)

**STEP 7. Variant Annotation**

Objective: Annotate the identified variants to determine their effect on genes and protein function.

Tool: [SnpEff](https://pcingola.github.io/SnpEff/) — a tool for annotating and predicting the effects of genetic variants.

We can also use [BCFtools `csq`](https://samtools.github.io/bcftools/bcftools.html#csq) for basic functional annotations.

*SnpEff Setup*:

##1. Download or build database for *E. coli* 
<pre>bash
java -jar snpEff.jar download Escherichia_coli_K12 </pre>

 ##2. Run annotation
<pre>bash
java -jar snpEff.jar Escherichia_coli_K12 raw_variants.vcf > annotated_variants.vcf </pre>

 Output: annotated_variants.vcf — annotated VCF with predicted effect

**Visualization Using IGV**

Objective: View aligned reads and variants on the reference genome to confirm variant accuracy and quality.

Tool: [IGV (Integrative Genomics Viewer)](https://software.broadinstitute.org/software/igv/) — a high-performance viewer for interactive exploration of genomics data.

Files to Load in IGV:
- `ecoli_reference.fasta` (reference genome)
- `aligned_sorted.bam` (sorted BAM file)
- `aligned_sorted.bam.bai` (index file)
- `annotated_variants.vcf` (final VCF file)


Steps:

##1. Launch IGV - Download and install from [https://software.broadinstitute.org/software/igv/download](https://software.broadinstitute.org/software/igv/download)

##2. Load Reference Genome
- Go to **Genomes > Load Genome from File**
- Select `ecoli_reference.fasta`

##3. Load Data Tracks
- Go to **File > Load from File**
- Load `aligned_sorted.bam` (IGV will auto-detect `.bai`)
- Load `annotated_variants.vcf`

##4. Navigate to Variants
- Use the **search bar** to jump to specific chromosome/position (e.g., `NC_000913.3:100000-100500`)
- Zoom in to examine read coverage and variant evidence

What we Looked For:
- Consistent read support for variants
- Read depth/coverage at variant sites
- Mapping quality and alignment confidence



**CONCLUSION**:
This visualization helps validate the bioinformatics pipeline by checking if variants are real (not artifacts) and well-supported by the data.

**KEY LEARNINGS**:
- Understanding of tools like **BWA**, **GATK**, **SAMtools**, **SnpEff**, and **IGV**
- Application of **variant calling techniques** to microbial genomics
- Gained insights into how computational pipelines can uncover **biological meaning from raw sequencing data**

>  Such variant analysis is critical in microbial research, outbreak tracing, epidemiology, and antibiotic resistance surveillance.

This project serves as a practical template for similar bacterial WGS analysis workflows in real-world research or diagnostics.
