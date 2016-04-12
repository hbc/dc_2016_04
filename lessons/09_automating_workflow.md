---
title: "Automating the workflow"
author: "Mary Piper, Meeta Mistry"
date: "Thursday, March 31, 2016"
---

Approximate time: 75 minutes

## Learning Objectives:

* Automate a workflow by grouping a series of sequential commands into a script
* Modify the script to grant it more flexibility 

## Calling Variants from all files

That was a lot of work, yes? But you have five more FASTQ files to go...

- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?


## Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...


### Collect all required commands
Using your command history, create a script file that will repeat these commands
for you. Name your script *run_variant_call_on_file.sh*. 

```bash

#! /bin/bash

### Alignmemt

# Align the reads using BWA
 bwa aln data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq/SRR097977.fastq_trim.fastq \
    > results/sai/SRR097977.aligned.sai

# Convert the .sai file to SAM
bwa samse data/ref_genome/ecoli_rel606.fasta results/sai/SRR097977.aligned.sai \
     data/trimmed_fastq/SRR097977.fastq > results/sam/SRR097977.aligned.sam

# Convert the SAM file to BAM
samtools view -S -b results/sam/SRR097977.aligned.sam >  results/bam/SRR097977.aligned.bam

# Sort the BAM file
samtools sort results/bam/SRR097977.aligned.bam results/bam/SRR097977.aligned.sorted.bam

### Variant calling

# Counting read coverage
samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
      results/bam/SRR097977.aligned.sorted.bam > results/bcf/SRR097977_raw.bcf

# Identify SNPs
bcftools view -bvcg results/bcf/SRR097977_raw.bcf > results/bcf/SRR097977_variants.bcf

# Filter SNPs
bcftools view results/bcf/SRR097977_variants.bcf \
       | /usr/share/samtools/vcfutils.pl varFilter - > results/vcf/SRR097977_final_variants.vcf

```

### Positional parameters (input files)
At the moment this script will only call variants for `SRR097977`, but we want to grant it more flexibility so that it can be run on any FASTQ file that we give it as input.

We can specify a filename as input using **positional parameters**. Positional parameters allow flexibility within a script.

The command-line arguments $1, $2, $3,...$9 are positional parameters, with $0 pointing to the actual command, program or shell script, and $1, $2, $3, ...$9 as the arguments to the command." This basically means that "Script Name" == $0, "First Parameter" == $1, "Second Parameter" == $2 and so on...

At the start of the script let's create a variable that captures an input parameter that must be supplied with the script name, when the user is running it. This input parameter will be the name of the file we want to work on. Let's also now specifiy where we will be getting the input data from:

```bash
# Get input file and locations  
fq=$1
data=data/trimmed_fastq
```

### Naming files (ouput files)
But, wait! What about the ouput files? They are all named using the SRA ID from the original FASTQ file. We can use a unix command to extract the base name of the file (without the path and .fastq extension) and ue this for naming all of the output files. We'll assign it to the `$base` variable:

```bash
# Grab base of filename for future naming
base=`basename $fq trim.fastq`
```

Now, wherever we see an output file that contains `SRR097977` we can substitute that with `$base`, and where the input file is required we replace that with `$fq`. For example:

```bash
# Align the reads using BWA
bwa aln data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq/$fq >  results/sai/$base.aligned.sai
```

### Assigning ouput file locations to variables
We can go a step further in flexibility. Since we've already created our output directories, we can now specify all of the paths for our output files in their proper locations. We will assign various **filenames and location to variables** both for convenience but also to make it easier to see what  is going on in the command below.

```bash

# set up output filenames and locations
sai=results/sai/$base.aligned.sai
sam=results/sam/$base.aligned.sam
bam=results/bam/${base}_aligned.bam
sorted_bam=results/bam/$base.aligned.sorted
raw_bcf=results/bcf/${base}_raw.bcf
variants=results/bcf/${base}_variants.bcf
final_variants=results/vcf/${base}_final_variants.vcf

```


### Genome reference variable
We can also add a **shortcut** to store the location to the **genome reference** FASTA file. Now for commands that use the genome (i.e bwa and samtools) we can run them using the genome variable (`$genome`) so these values are not static and hard-coded.

```bash
# location to genome reference FASTA file
genome=data/ref_genome/ecoli_rel606.fasta
```


### Usage and commenting
Finally, when you are writing a multi-step workflow that accepts command-line options (positional parameters), it is very important to **have the usage right in the beginning of the script**. In addition to the usage, it is a good practice to comment about the inputs, steps/tools and output briefly. It is easier to fill both of these in after your script is ready and you have done a test run. So even though this is at the top of the script, it may be the last few lines you add to the script.

Now, walking thru the code, we can make some substitutions as described above. We will also add a command at the beginning of the script to change directories into `~/dc_workshop` so that this script can be run from any location.


### Final script

```bash

#!/bin/bash

# USAGE: sh run_variant_call_on_file.sh <fastq file> 
# This script will take the location and name of a fastq file and perform the following steps on it in a new directory. 
    ## starting with alignment using bwa-backtrack
    ## convert the .sai to SAM -> BAM -> sorted BAM, 
    ## counting read coverage with samtools, 
    ## identify SNPs and filter SNPs.


cd ~/dc_workshop

# Get input file and locations  
fq=$1
data=data/trimmed_fastq

# Grab base of filename for future naming
base=`basename $fq .fastq_trim.fastq`
echo "Starting analysis of" $base

# set up output filenames and locations
sai=results/sai/$base.aligned.sai
sam=results/sam/$base.aligned.sam
bam=results/bam/${base}_aligned.bam
sorted_bam=results/bam/$base.aligned.sorted
raw_bcf=results/bcf/$base_raw.bcf
variants=results/bcf/$base_variants.bcf
final_variants=results/vcf/$base_final_variants.vcf

# location to genome reference FASTA file
genome=data/ref_genome/ecoli_rel606.fasta

### Alignmemt

# Align the reads using BWA
bwa aln $genome $data/${base}.fastq_trim.fastq > $sai

# Convert the .sai file to SAM
bwa samse $genome $sai $data/${base}.fastq_trim.fastq > $sam

# Convert the SAM file to BAM
samtools view -S -b $sam >  $bam

# Sort the BAM file
samtools sort $bam $sorted_bam

### Variant calling

# Counting read coverage
samtools mpileup -g -f $genome ${sorted_bam}.bam > $raw_bcf

# Identify SNPs
bcftools view -bvcg $raw_bcf > $variants

# Filter SNPs
bcftools view $variants | /usr/share/samtools/vcfutils.pl varFilter - > $final_variants
```








