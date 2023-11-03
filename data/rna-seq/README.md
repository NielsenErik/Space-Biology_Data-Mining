# RNA-seq data for space project

### Columns

In the file GLDS-48_rna_seq_Normalized_Counts.csv, each column
represents one sample, in our case, one mouse is one sample. We are
interested in following mice, all other columns in the file can be
ignored.

  -----------------------------------------------------------------------
  Name of the mouse                   Environment
  ----------------------------------- -----------------------------------
  Mmus_C57-6J_LVR_GC_C_Rep1_M36       Ground Control

  Mmus_C57-6J_LVR_GC_C_Rep2_M37       Ground control

  Mmus_C57-6J_LVR_GC_C_Rep3_M38       Ground Control

  Mmus_C57-6J_LVR_GC_C_Rep4_M39       Ground control

  Mmus_C57-6J_LVR_GC_C_Rep5_M40       Ground Control

  Mmus_C57-6J_LVR_FLT_C_Rep1_M25      Space flight

  Mmus_C57-6J_LVR_FLT_C_Rep2_M26      Space flight

  Mmus_C57-6J_LVR_FLT_C_Rep3_M27      Space flight

  Mmus_C57-6J_LVR_FLT_C_Rep4_M28      Space flight

  Mmus_C57-6J_LVR_FLT_C_Rep5_M30      Space flight
  -----------------------------------------------------------------------

### Rows

Each row in the file represents one gene. The first column contains the
Ensembl gene IDs for these genes. If you want to know more about a gene.
Search for the Ensembl gene ID in the Ensembl database:
[[https://www.ensembl.org/index.html]{.underline}](https://www.ensembl.org/index.html)

For example:

Using the Ensembl database I learn that the gene ENSMUSG00000000441
which is present in our file is also called Raf1, or v-raf-leukemia
viral oncogene 1.

### What does the numbers in the file actually mean?

This file contains normalised counts. What does that mean?

After you have performed an RNA-seq experiment (added the samples to the
complicated machine in the lab) you get a computer file with all the RNA
that was present in the sample.

Unfortunately, the universe has not annotated the RNA, so we don\'t know
which is which. In order to know the names of the genes which have
encoded these RNA, the data needs to be *aligned*.

After the data has been aligned, you get a file with the *raw counts*. A
raw counts table gives us how many RNA was present for each gene in a
sample.

So lets say that in sample 1 we find 10 RNA from gene A and in sample 2
we find 5 RNA from gene A. This should mean that in sample 1 we have a
higher activation of gene A than in sample 2 right?

Not necessarily, sample 2 might have been more diluted. To compare
RNA-seq data from different samples you need to normalise the data.
Meaning that instead of just counting the number of RNA, you present it
as a percentage: **"In sample 1 we found that gene A has 7 RNA out of
100 million total RNA."**

The data we are using was normalised using the tool DESeq2. More about
DESeq2 can be found here:
[[https://chipster.csc.fi/manual/deseq2.html]{.underline}](https://chipster.csc.fi/manual/deseq2.html)

I think that the data presented in our file is in the unit *RPKM* (reads
per kilobase of exon per million reads mapped). What does that mean? It
doesn\'t matter. High number means that the gene is activated a lot, low
number means it\'s less activated.

### So what are we using for our analysis?

Since we have multiple samples (5 from each environment), we might just
want to add all the samples into one and take the average for each gene,
but I will do some research regarding if this actually is the best idea.

Maybe we need to visualise the data somehow in order to find outliers.

**Protocol used for RNA-seq data processing in the researches own
words**

Raw fastq files were filtered using Trim Galore! (version 0.6.2).
Trimmed fastq file quality was evaluated with FastQC (version 0.11.8),
and MultiQC (version 1.7) was used to generate MultiQC reports. Mus
musculus STAR and RSEM references were built using STAR (version 2.7.1a)
and RSEM (version 1.3.1), respectively, Ensembl release 96, genome
version mm10-GRCm38 (Mus_musculus.GRCm38.dna.toplevel.fa), and the
following gtf annotation file: Mus_musculus.GRCm38.96.gtf. Trimmed reads
were aligned to the Mus musculus STAR reference with STAR (version
2.7.1a) and aligned reads were quantified using RSEM (version 1.3.1).
The runsheet was generated with dp_tools (version 1.1.8) and the
runsheet and quantification data were imported to R (version 4.1.3) with
tximport (version 1.27.1) and normalized with DESeq2 (version 1.34.0).
Differential expression analysis was performed in R (version 4.1.3)
using DESeq2 (version 1.34.0); all groups were compared using the Wald
test and the likelihood ratio test was used to generate the F statistic
p-value. Gene annotations were assigned using the GeneLab Reference
Annotations pipeline
(https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110.md),
with STRINGdb (version 2.8.4), PANTHER.db (version 1.0.11), and
org.Mm.eg.db (version 3.15.0). Processing scripts used can be found on
GitHub here:
https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/GLDS_Processing_Scripts/GLDS-48/GLDS_version_11

Document written by Zacharina 2023-11-02

Last updated: 2023-11-02