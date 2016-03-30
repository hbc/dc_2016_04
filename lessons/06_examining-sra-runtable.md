---
layout: lesson
root: .
title: Examining data on the NCBI SRA database
minutes: 20
---

## Learning Objectives

* Get a general understanding of the NCBI SRA
* Understand how to download a summary run table of the SRA
* Be able to open a run table using a spreadsheet editor

## Lesson

The NCBI Sequence Read Archive is a large (>3 quadrillion basepairs as of 2014) repository for next-generation sequence data. Like many NCBI database it is complex and mastering its use is greater than the scope of this lesson. Very often, as in the Lenski paper, there will be a direct link (perhaps in the supplemental information) to where on the SRA the dataset can be found. The link from the Lenski paper is: [http://www.ncbi.nlm.nih.gov/sra?term=SRA026813](http://www.ncbi.nlm.nih.gov/sra?term=SRA026813)



This section takes place on your own laptop - no need to open the shell or connect to your remote computer yet!

#### A. Locate the Run Accessor for the Lenski Dataset on the SRA

1. Access the Lenski dataset from the provided link: [http://www.ncbi.nlm.nih.gov/sra?term=SRA026813](http://www.ncbi.nlm.nih.gov/sra?term=SRA026813).  
You will be presented with a page for the overall SRA accession SRA026813 - this is a collection of all the experimental data
2. Click on the first entry ([ZDB30](http://www.ncbi.nlm.nih.gov/sra/SRX040669%5Baccn%5D)); this will take you to a page for an SRX (Sequence Read eXperiment). Take a few minutes to examine some of the descriptions on the page
3. Click on the ['All runs'](http://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP004752) link under where it says **Study**; this is a description of all of the NGS datasets related to the experiment.
4. Go to the top of the page and in the **Total** row you will see there are 37 runs, 10.15Gb data, and 16.45 Gbases of data. Click the 'RunInfo Table' button.

We are not downloading any actual sequence data here! This is only a text file that fully describes the entire dataset

You should now have a file called **SraRunTable.txt**

#### B. Review the SraRunTable in a spreadsheet program


1. Using your choice of spreadsheet program open the **SraRunTable.txt** file. If prompted this is a tab-delimited file.

### Exercise

Discuss with the person next to you:

1. What strain of *E. coli* was used in this experiment?
2. What was the sequencing platform used for this experiment?
3. What samples in the experiment contain [paired end](http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html) sequencing data?
4. What other kind of data is available?
5. Are you collecting this kind of information about your sequencing runs?

After answering the question, you should avoid saving this file; we don't want to make any changes. If you were to save this file, make sure you save it as a plain **.txt** file.


## Where to learn more

#### About the Sequence Read Archive

* You can learn more about the SRA by reading the [SRA Documentation](http://www.ncbi.nlm.nih.gov/Traces/sra/)
* The best way to transfer a large SRA dataset is by using the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc)

Previous: [Introducing the dataset](01-intro-to-dataset.html)
