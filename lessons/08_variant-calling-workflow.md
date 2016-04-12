---
title: "Variant calling workflow"
author: "Mary Piper, Meeta Mistry"
date: "Thursday, March 31, 2016"
---

Approximate time: 120 minutes

## Learning Objectives:

* Understand the different steps involved in variant calling
* Use a series of command line tools to execute a variant calling workflow
* Becoming familiar with data formats encountered during variant calling


## Setting up

To get started with this lesson, make sure you are in `dc_workshop`. Now let's copy over the reference data required for alignment:

```bash
$ cd ~/dc_workshop
$ cp -r ~/.dc_sampledata_lite/ref_genome/ data/
```

Your directory structure should now look like this:

<pre>
dc_workshop
├── data
    ├── ref_genome
        └── ecoli_rel606.fasta
    ├── untrimmed_fastq
    └── trimmed_fastq
        ├── SRR097977.fastq
        ├── SRR098026.fastq
        ├── SRR098027.fastq
        ├── SRR098028.fastq
        ├── SRR098281.fastq
        └── SRR098283.fastq
 ├── results
 └── docs

</pre>

You will also need to create directories for the results that will be generated as part of the workflow: 
```bash
$ mkdir  results/sai results/sam results/bam results/bcf results/vcf
```

> *NOTE: All of the tools that we will be using in this workflow have been pre-installed on our remote computer*

## Alignment to a reference genome

We have already trimmed our reads so now the next step is alignment of our quality reads to the reference genome.

![workflow_align](../img/variant_calling_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to choose from and while there is no gold standard there are some tools that are better suited for particular NGS analyses. We will be using the [Burrows Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent sequences against a large reference genome. The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

### Index the reference genome
Our first step is to index the reference genome for use by BWA. *NOTE: This only has to be run once*. The only reason you would want to create a new index is if you are working with a different reference  genome or you are using a different tool for alignment.

```bash    
$ bwa index data/ref_genome/ecoli_rel606.fasta     # This step helps with the speed of alignment
```

Eventually we will loop over all of our files to run this workflow on all of our samples, but for now we're going to work on just one sample in our dataset `SRR098283.fastq`:

```bash
$ ls -alh ~/dc_workshop/data/trimmed_fastq/SRR098283.fastq_trim.fastq 
```

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an aligner. BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate.

Since we are working with short reads we will be using BWA-backtrack. The usage for BWA-backtrack is `bwa aln path/to/ref_genome.fasta path/to/fastq > SAIfile`. This will create a `.sai` file which is an intermediate file containing the suffix array indexes. 
    
Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default parameters here, your use case might require a change of parameters. *NOTE: Always read the manual page for any tool before using and try to understand the options.*

```bash
$ bwa aln data/ref_genome/ecoli_rel606.fasta \
    data/trimmed_fastq/SRR098283.fastq_trim.fastq > results/sai/SRR098283.aligned.sai
```

## Alignment cleanup

![workflow_clean](../img/variant_calling_workflow_cleanup.png)

Post-alignment processing of the alignment file includes:

1. Converting output SAI alignment file to a BAM file
2. Sorting the BAM file by coordinate

### Convert the format of the alignment to SAM/BAM

The SAI file is not a standard alignment output file and will need to be converted into a SAM file before we can do any downstream processing. 

#### SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam)), is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not have time to go in detail of the features of the SAM format, the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification. **The binary version of SAM is called a BAM file.**

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.

![sam_bam1](../img/sam_bam.png)


![sam_bam2](../img/sam_bam3.png)

First we will use the `bwa samse` command to convert the .sai file to SAM format:

```bash
$ bwa samse data/ref_genome/ecoli_rel606.fasta \
      results/sai/SRR098283.aligned.sai \
      data/trimmed_fastq/SRR098283.fastq_trim.fastq > \
      results/sam/SRR098283.aligned.sam
```
Explore the information within your SAM file:

```bash
$ head results/sam/SRR098283.aligned.sam
```	
Now convert the SAM file to BAM format for use by downstream tools: 

```bash
$ samtools view -S -b results/sam/SRR098283.aligned.sam > results/bam/SRR098283.aligned.bam
```
### Sort BAM file by coordinates

Sort the BAM file:

```bash
$ samtools sort -f results/bam/SRR098283.aligned.bam results/bam/SRR098283.aligned.sorted.bam
```

*SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.*


## Variant calling

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given position in an individual genome or transcriptome, often referred to as a Single Nucleotide Polymorphism (SNP). The call is usually accompanied by an estimate of variant frequency and some measure of confidence. Similar to other steps in this workflow, there are number of tools available for variant calling. In this workshop we will be using `bcftools`, but there are a few things we need to do before actually calling the variants.

![workflow](../img/variant_calling_workflow.png)

### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with samtools [mpileup](http://samtools.sourceforge.net/mpileup.shtml):

    $ samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
    results/bam/SRR098283.aligned.sorted.bam > results/bcf/SRR098283_raw.bcf

***We have only generated a file with coverage information for every base with the above command; to actually identify variants, we have to use a different tool from the samtools suite called [bcftools](https://samtools.github.io/bcftools/bcftools.html).***

### Step 2: Detect the single nucleotide polymorphisms (SNPs)

Identify SNPs using bcftools:

    $ bcftools view -bvcg results/bcf/SRR098283_raw.bcf > results/bcf/SRR098283_variants.bcf

### Step 3: Filter and report the SNP variants in VCF (variant calling format)

Filter the SNPs for the final output in VCF format, using vcfutils.pl:

    $ bcftools view results/bcf/SRR098283_variants.bcf | /usr/share/samtools/vcfutils.pl varFilter - > 		results/vcf/SRR098283_final_variants.vcf
	
*`bcftools view` converts the binary format of bcf files into human readable format (tab-delimited) for `vcfutils.pl` to perform the filtering. Note that the output is in VCF format, which is a text format.*

## Explore the VCF format:

	$ less results/vcf/SRR098283_final_variants.vcf

You will see the **header** which describes the format, when the file was created, the tools version along with the command line parameters used and some additional column information:

	##reference=file://data/ref_genome/ecoli_rel606.fasta
	##contig=<ID=NC_012967.1,length=4629812>
	##ALT=<ID=X,Description="Represents allele(s) other than observed.">
	##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
	##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
	##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
	.
	.
	.
	.
	##bcftools_callVersion=1.2+htslib-1.2.1
	##bcftools_callCommand=call -cv -O b results/bcf/SRR098283_raw.bcf
	##bcftools_viewVersion=1.2+htslib-1.2.1
	##bcftools_viewCommand=view results/bcf/SRR098283_variants.bcf

Followed by the **variant information**:

	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  results/bam/SRR098283.trimmed.aligned.sorted.bam
	NC_012967.1     110152  .       T       A       18.0963 .       DP=3;VDB=0.74;SGB=-0.453602;RPB=1;MQB=1;MQSB=1;BQB=1;MQ0F=0;AF1=0.502509;AC1=1;DP4=1,0,0,2;MQ=37;FQ=-7.78372;PV4=0.333333,1,1,0.20326   GT:PL   0/1:48,0,20
	NC_012967.1     270633  .       G       T       26.7735 .       DP=2;VDB=0.76;SGB=-0.453602;MQ0F=0;AF1=1;AC1=2;DP4=0,0,2,0;MQ=37;FQ=-32.988     GT:PL   1/1:58,6,0
	NC_012967.1     475173  .       G       C       21.7931 .       DP=2;VDB=0.14;SGB=-0.453602;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,1,1;MQ=37;FQ=-32.988      GT:PL   1/1:53,6,0
	NC_012967.1     1017485 .       G       T       5.46014 .       DP=5;VDB=0.58;SGB=-0.453602;RPB=0;MQB=1;MQSB=1;BQB=0.5;MQ0F=0;AF1=0.499901;AC1=1;DP4=0,3,2,0;MQ=37;FQ=7.77964;PV4=0.1,0.268358,1,1      GT:PL   0/1:34,0,56

The first columns represent the information we have about a predicted variation. 

CHROM and POS provide the config information and position where the variation occurs. 

ID is a `.` until we add annotation information. 

REF and ALT represent the genotype at the reference and in the sample, always on the foward strand. 

QUAL then is the Phred scaled probablity that the observed variant exists at this site. Ideally you would need nothing else to filter out bad variant calls, but in reality we still need to filter on multiple other metrics. 

The FILTER field is a `.`, i.e. no filter has been applied, otherwise it will be set to either PASS or show the (quality) filters this variant failed. 



The last columns contains the genotypes and can be a bit more tricky to decode. In brief, we have:

* GT: The genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc.

* GQ: the Phred-scaled confidence for the genotype

* AD, DP: Reflect the depth per allele by sample and coverage

* PL: the likelihoods of the given genotypes

The BROAD's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place to learn more about VCF file format.


## Assess the alignment (visualization) - optional step

Index the BAM file for visualization with IGV:

    $ samtools index results/bam/SRR098283.aligned.sorted.bam

**Transfer files to your laptop**

Using FileZilla, transfer the following 3 files to your local machine, 
`results/bam/SRR098283.aligned.sorted.bam`,

`results/bam/SRR098283.aligned.sorted.bam.bai`, 

`data/ref_genome/ecoli_rel606.fasta`

**Visualize**	

* Start [IGV](https://www.broadinstitute.org/software/igv/download)
* Load the genome file into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*


