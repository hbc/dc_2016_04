---
layout: lesson
root: .
title: Introduction to workshop dataset
minutes: 45
---

## Learning Objectives
* Understand how and why we choose this dataset
* Have a general idea of the experiment and its objectives
* Get a general understanding of the NCBI SRA
* Understand how to download a summary run table of the SRA
* Be able to open a run table using a spreadsheet editor

## Introduction to the dataset

#### Description of dataset

Microbes are ideal organisms for exploring 'Long-term Evolution Experiments' (LTEEs) - thousands of generations can be generated and stored in a way that would be virtually impossible for more complex eukaryotic systems. In [Blount et al 2012](Lenski_paper.pdf), 12 populations of *Escherichia coli* were propagated for more than 40,000 generations in a glucose-limited minimal medium. This medium was supplemented with citrate which *E. coli* cannot metabolize in the aerobic conditions of the experiment. Sequencing of the populations at regular time points reveals that spontaneous citrate-using mutants (Cit+) appeared in a population of *E.coli* (designated Ara-3) at around 31,000 generations. It should be noted that spontaneous Cit+ mutants are extraordinarily rare - inability to metabolize citrate is one of the defining characters of the *E. coli* species. Eventually, Cit+ mutants became the dominant population as the experimental growth medium contained a high concentration of citrate relative to glucose.

Strains from generation 0 to generation 40,000 were sequenced, including ones that were both Cit+ and Cit- after generation 31,000.

For the purposes of this workshop we're going to be working with 6 of the samples from this experiment. We also made up genome sizes for each of the strains, to look at the relationship between Cit status and genome size.  **The genome sizes are not real data!!**


| SRA Run Number | Clone | Generation | Cit | GenomeSize |
| -------------- | ----- | ---------- | ----- | ----- |
| SRR098028 | REL1166A | 2,000 | Unknown | 4.63 |
| SRR098281 | ZDB409 | 5,000 | Unknown | 4.6 |
| SRR098283 | ZDB446 | 15,000 | Cit- | 4.66 |
| SRR097977 | CZB152 | 33,000 | Cit+ | 4.8 |
| SRR098026 | CZB154 | 33,000 | Cit+ | 4.76 |
| SRR098027 | CZB199 | 33,000 | Cit- | 4.59 |


We want to be able to look at the genome size to see if there is a difference between genome size and the Cit status of the strain. We also want to analyze the sequences to figure out what changes occurred in genomes to make the strains Cit+. Ultimately, we will answer the questions:

- What is the distribution of genome sizes for all the strains?
- Is there a relationship between genome size and Cit status?
- How many base pair changes are there between the Cit+ and Cit- strains?
- What are the base pair changes between strains?

#### Accessing the original archived data
The sequencing dataset was attained from the [NCBI Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra), which is a large (>3 quadrillion basepairs as of 2014) repository for next-generation sequence data. Like many NCBI databases, it is complex and mastering its use is greater than the scope of this lesson. Very often, as in the Lenski paper, there will be a direct link (perhaps in the supplemental information) to where on the SRA the dataset can be found. The link from the Lenski paper is: [http://www.ncbi.nlm.nih.gov/sra?term=SRA026813](http://www.ncbi.nlm.nih.gov/sra?term=SRA026813)

###### A. Locate the Run Accessory Table for the Lenski Dataset on the SRA

1. Access the Lenski dataset from the provided link: [http://www.ncbi.nlm.nih.gov/sra?term=SRA026813](http://www.ncbi.nlm.nih.gov/sra?term=SRA026813).  
You will be presented with a page for the overall SRA accession SRA026813 - this is a collection of all the experimental data
2. Click on the first entry ([ZDB30](http://www.ncbi.nlm.nih.gov/sra/SRX040669%5Baccn%5D)); this will take you to a page for an SRX (Sequence Read eXperiment). Take a few minutes to examine some of the descriptions on the page
3. Click on the ['All runs'](http://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP004752) link under where it says **Study**; this is a description of all of the NGS datasets related to the experiment.
4. Go to the top of the page and in the **Total** row you will see there are 37 runs, 10.15Gb data, and 16.45 Gbases of data. Click the 'RunInfo Table' button.

We are not downloading any actual sequence data here! This is only a text file that fully describes the entire dataset

You should now have a file called **SraRunTable.txt**

###### B. Review the SraRunTable in a spreadsheet program


1. Using your choice of spreadsheet program open the **SraRunTable.txt** file. If prompted this is a tab-delimited file.

***
**Exercises**

Discuss with the person next to you:

1. What strain of *E. coli* was used in this experiment?
2. What was the sequencing platform used for this experiment?
3. What samples in the experiment contain [paired end](http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html) sequencing data?
4. What other kind of data is available?
5. Are you collecting this kind of information about your sequencing runs?
***
After answering the questions, you should avoid saving this file; we don't want to make any changes. If you were to save this file, make sure you save it as a plain **.txt** file.


## Where to learn more

#### About the Sequence Read Archive

* You can learn more about the SRA by reading the [SRA Documentation](http://www.ncbi.nlm.nih.gov/Traces/sra/)
* The best way to transfer a large SRA dataset is by using the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc)

#### References

Blount, Z.D., Barrick, J.E., Davidson, C.J., Lenski, R.E.
Genomic analysis of a key innovation in an experimental Escherichia coli population (2012) Nature, 489 (7417), pp. 513-518.  
[Paper](Lenski_paper.pdf), [Supplemental materials](Lenski-s1.pdf)  
Data on NCBI SRA: [http://www.ncbi.nlm.nih.gov/sra?term=SRA026813](http://www.ncbi.nlm.nih.gov/sra?term=SRA026813)


