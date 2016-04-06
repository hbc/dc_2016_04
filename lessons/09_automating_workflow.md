
---
title: "Automating the workflow"
author: "Mary Piper, Meeta Mistry"
date: "Thursday, March 31, 2016"
---

Approximate time: 75 minutes

## Learning Objectives:

* Automate a workflow by grouping a series of sequential commands into a script
* Modify the script to grant it more flexibility 

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
