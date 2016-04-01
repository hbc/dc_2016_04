# Lesson

Automating a workflow
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?

* Use a series of command line tools to execute a variant calling workflow
* Automate a workflow by grouping a series of sequential commands into a script
* Modify and submit the workflow script to the cluster

# Running a Workflow

## Setting up

To get started with this lesson, we will need to grab some data from an outside
server using `wget` on the command line.

Make sure you are in the dc_workshop directory first

```bash
$ cd ~/dc_workshop
$ wget http://reactomerelease.oicr.on.ca/download/archive/variant_calling.tar.gz
```

The file 'variant_calling.tar.gz' is what is commonly called a "tarball", which is
a compressed archive similar to the .zip files we have seen before.  We can decompress
this archive using the command below.

```bash
$ tar -zxvf variant_calling.tar.gz
```
This will create a directory tree that contains some input data (reference genome and fastq files)
and a shell script that details the series of commands used to run the variant calling workflow.

<pre>
variant_calling
├── data
    ├── ref_genome
        └── ecoli_rel606.fasta
    └── trimmed_fastq
        ├── SRR097977.fastq
        ├── SRR098026.fastq
        ├── SRR098027.fastq
        ├── SRR098028.fastq
        ├── SRR098281.fastq
        └── SRR098283.fastq
└── run_variant_calling.sh
</pre>

Remember the variant calling workflow.

![workflow](../img/variant_calling_workflow.png)

We'll first perform the commands for all the above steps (run through the workflow). 

Next, we'll create a script for the commands and test it. 

Finally, we'll modify the script to run on the cluster.

So let's get started.

Change directories to the `variant_calling` directory:

	$ cd ~/dc_workshop/variant_calling     

Create the directory structure for our results:

	$ mkdir -p results/sai results/sam results/bam results/bcf results/vcf

## Alignment to a reference genome

Since the quality control steps have already been completed, and we are starting with trimmed reads, we can start our workflow with alignment of our quality reads to the reference genome.

![workflow_align](../img/variant_calling_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

### Index the reference genome
Our first step is to index the reference genome for use by BWA. BWA
is a program that is pre-installed on our server:
    
	bwa index data/ref_genome/ecoli_rel606.fasta     # This step helps with the speed of alignment
	
In the script, we will eventually loop over all of our files and have the cluster work
on each one in parallel. For now, we're going to work on just one to set up our workflow:

```bash
ls -alh ~/dc_workshop/variant_calling/data/trimmed_fastq/SRR098283.fastq
```

### Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against, and performing the read alignment using one of several alignment tools such as [NovoAlign](http://www.novocraft.com/main/page.php?s=novoalign) or [BWA](https://github.com/lh3/bwa). 

The usage for BWA is `bwa aln path/to/ref_genome.fasta path/to/fastq > SAIfile`
    
Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). While we are running bwa with the default parameters here, your use case might require a change of parameters. NOTE: Always read the manual page for any tool before using and try to understand the options.

    bwa aln data/ref_genome/ecoli_rel606.fasta \
    data/trimmed_fastq/SRR098283.fastq > \
    results/sai/SRR098283.aligned.sai

## Alignment cleanup

![workflow_clean](../img/variant_calling_workflow_cleanup.png)

Post-alignment processing of the alignment file includes:

1. Converting output SAI alignment file to a BAM file
2. Sorting the BAM file by coordinate

### Convert the format of the alignment to BAM

Before we convert to BAM format, we need to first convert the output to the SAM format. Usage: `bwa samse path/to/ref_genome.fasta path/to/SAIfile path/to/fastq > SAMfile`

	bwa samse data/ref_genome/ecoli_rel606.fasta \
	results/sai/SRR098283.aligned.sai \
	data/trimmed_fastq/SRR098283.fastq > \
	results/sam/SRR098283.aligned.sam

Explore the information within a [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam):

	head results/sam/SRR098283.aligned.sam
	
Now convert the SAM file to BAM format for use by downstream tools: 

    samtools view -S -b results/sam/SRR098283.aligned.sam > \
    results/bam/SRR098283.aligned.bam

### Sort BAM file by coordinates -- doesn't work!

Sort the BAM file:

    samtools sort -O bam -T temp.prefix results/bam/SRR098283.aligned.bam > \
    results/bam/SRR098283.aligned.sorted.bam

*SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.*

### Assess the alignment (visualization) - optional step

Index the BAM file for visualization with IGV:

    samtools index results/bam/SRR098283.aligned.sorted.bam

**Transfer files to your laptop**

Using FileZilla, transfer the following 3 files to your local machine, 
`results/bam/SRR098283.aligned.sorted.bam`,

`results/bam/SRR098283.aligned.sorted.bam.bai`, 

`data/ref_genome/ecoli_rel606.fasta`

**Visualize**	

* Start [IGV](https://www.broadinstitute.org/software/igv/download)
* Load the genome file into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*

## Call variants

![workflow](../img/variant_calling_workflow.png)

### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with samtools [mpileup](http://samtools.sourceforge.net/mpileup.shtml):

	samtools faidx data/ref_genome/ecoli_rel606.fasta     # We will a reference genome index using samtools for variant calling as well 

    samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
      results/bam/SRR098283.aligned.sorted.bam > results/bcf/SRR098283_raw.bcf

***We have only generated a file with coverage information for every base with the above command; to actually identify variants, we have to use a different tool from the samtools suite called [bcftools](https://samtools.github.io/bcftools/bcftools.html).***

### Step 2: Detect the single nucleotide polymorphisms (SNPs)

Identify SNPs using bcftools:

    bcftools call -vc -O b results/bcf/SRR098283_raw.bcf > results/bcf/SRR098283_variants.bcf

### Step 3: Filter and report the SNP variants in VCF (variant calling format)

Filter the SNPs for the final output in VCF format, using vcfutils.pl:

    bcftools view results/bcf/SRR098283_variants.bcf | vcfutils.pl varFilter - > \
    results/vcf/SRR098283_final_variants.vcf
	
*`bcftools view` converts the binary format of bcf files into human readable format (tab-delimited) for `vcfutils.pl` to perform the filtering. Note that the output is in VCF format, which is a text format.*

Explore the VCF format:

	less results/vcf/SRR098283_final_variants.vcf

You will see the header which describes the format, when the file was created, the tools version along with the command line parameters used and some additional column information:

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

Followed by the variant information:

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

#### Exercise I - Calling Variants from all files?

That's a lot of work, yes? But you have five more FASTQ files to go...

- Try running this workflow on a different FASTQ file. What did you have to do differently
in order to get this workflow to work?
- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?


#### Exercise II - Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...

- Using your command history, create a script file that will repeat these commands
for you. Name your script *run_variant_call_on_file.sh*. Delete your results 
directories, and run your script. Do you get all the proper output files?

One additional command we can put in the top of the script to allow you to see what
is going on is the `set -x` bash command. This debugging tool will print out every
step before it is executed.

- Insert the debugging command in your script and re-run it. How is the output different?
If you're comfortable with how this looks and runs, then comment out this line.

- In order run this workflow on another file, you'll need to make changes. Copy this file,
giving the file a similar name, and make appropriate changes to run on another input
FASTQ file. What did you have to do differently in order to get this workflow to work?

- Knowing techniques that you've learned in previous lessons, what can we do to make this
workflow more friendly to different sets of input files?

- Again, reviewing your two scripts, are there additional commonalities across scripts
or within scripts that we could optimize?


#### Exercise III - Granting our Workflow More Flexibility

A couple of changes need to be made to make this script more friendly to both changes
in the workflow and changes in files. 

The first major change is allowing a change in the filename. Thus at the start of 
the script let's capture an input parameter that must be supplied with the script name.
This input parameter will be the name of the file we want to work on:

     fq="$1"

And we'll add a shortcut to store the location to the genome reference FASTA file:

     # location to genome reference FASTA file
     genome=data/ref_genome/ecoli_rel606.fasta


Now, walking thru the code, we can make some substitutions. To index with bwa and samtools
we can run those commands on the genome variable ($genome) so these values aren't 
static & hardcoded:

     bwa index $genome
     samtools faidx $genome

We'll keep the output paths creation, as it looks fine. (Though really, we could
put results/ in a variable and declare that at the top, so we can change where the
results will be as well. We'll leave that for an optional exercise)

     # make all of our output directories
     mkdir -p results/sai
     mkdir -p results/sam
     mkdir -p results/bam
     mkdir -p results/bcf
     mkdir -p results/vcf

In the script, it is a good idea to use `echo` for debugging/reporting to the screen

    echo "Processing file $fq ..."

We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it
to the $base variable

    # grab base of filename for future naming
    base=$(basename $fq .fastq)
    echo "basename is $base"

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the sommand below.

    # set up output filenames and locations
    fq=data/trimmed_fastq/$base\.fastq
    sai=results/sai/$base\_aligned.sai
    sam=results/sam/$base\_aligned.sam
    bam=results/bam/$base\_aligned.bam
    sorted_bam=results/bam/$base\_aligned_sorted.bam
    raw_bcf=results/bcf/$base\_raw.bcf
    variants=results/bcf/$base\_variants.bcf
    final_variants=results/vcf/$base\_final_variants.vcf

Our data are now staged.  We now need to change the series of command below
to use our variables so that it will run with more flexibility the steps of the 
analytical workflow:

    # Align the reads to the reference genome
    bwa aln $genome $fq > $sai

    # Convert the output to the SAM format
    bwa samse $genome $sai $fq > $sam

    # Convert the SAM file to BAM format
    samtools view -S -b $sam > $bam

    # Sort the BAM file
    samtools sort -O 'bam' -T temp.prefix $bam $sorted_bam

    # Index the BAM file for display purposes
    samtools index $sorted_bam

    # Do the first pass on variant calling by counting read coverage
    samtools mpileup -g -f $genome $sorted_bam > $raw_bcf

    # Do the SNP calling with bcftools
    bcftools call -vc -O b $raw_bcf > $variants

    # And finally, filter the SNPs for the final output
    bcftools view $variants | vcfutils.pl varFilter - > $final_variants

This new script is now ready for running:
	
	sh run_variant_call_on_file.sh <name of fastq>

#### Exercise IV - Parallelizing workflow for efficiency - NEEDS TO CHANGE

To run the same script on a worker node on the cluster via the job scheduler, we need to add our **SLURM directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the
compute node(s). 

Copy the `run_variant_call_on_file.sh` file and give it a new name `run_variant_call_on_file.sbatch`. Add the SLURM directives to the beginning of the file.

So the top of the file should look like:

    #!/bin/bash
    #
    #SBATCH -p serial_requeue   # Partition to submit to (comma separated)
    #SBATCH -n 1                # Number of cores
    #SBATCH -N 1                # Ensure that all cores are on one machine
    #SBATCH -t 0-1:00           # Runtime in D-HH:MM (or use minutes)
    #SBATCH --mem 100           # Memory in MB
    #SBATCH -J var_call_ecoli      # Job name
    #SBATCH -o var_call_ecoli.out       # File to which standard out will be written
    #SBATCH -e var_call_ecoli.err       # File to which standard err will be written
    #SBATCH --mail-type=ALL     # Type of email notification: BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=<your-email@here.com> # Email to which notifications will be sent 

What we'd like to do is run this script on a compute node for every trimmed FASTQ -- pleasantly parallelizing our workflow. And now is where we'll use the for loop with the power of the cluster: 

    for fq in data/trimmed_fastq/*.fastq
    do
      sbatch run_variant_call_on_file.sh $fq
      sleep 1
    done

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that you submitted.

You can see their progress by using the squeue command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the scancel command, should something go wrong and you need to
cancel your jobs.

* Change the script so that one can include an additional variable to point to
 a results directory.
