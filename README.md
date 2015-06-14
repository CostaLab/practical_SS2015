# costalab
Bioinformatics Lab RWTH Aachen SS '15

# Running alignments

The script bio_pipeline.sh allows an easy and almost-painless alignment procedure.

At the very least, it requires two files:

 * a FASTA reference genome (option -ref) (this MUST be in the current directory)
 * a SRA reads file (option -sra) (this can be in any directory, or even "fake" with option -sra-nocheck)

Without other options, the script will:

 * assume the organism is haploid (ie, only has one set of chromosomes)
 * extract the SRA file to FASTQ files, and remove (in simple cases) the adaptor-only files
  * the previous step also infers whether the SRA is paired-ends or single-ends
 * align the reads to the reference genome, using bwa-backtrack with default options
 * convert SAM file to BAM
 * remove duplicates, re-align near indels, produce a clean BAM file
 * generate plenty of statistics for the final BAM file
 * run GATK SNP calling and print the final number of SNPs

All the files are kept: it is up to you to delete the temporary files, or the unnecessary ones. A log is generated to output.log, containing both the stdout and stderr produced during the execution of bio_pipeline.sh.

Among the options, some are notable:

 * -mem: use bwa-mem instead of backtrack (very important for Ion Torrent, 454)
 * -mem-pacbio: a special option for PacBio, instead of the previous one (adds "-x pacbio" to the bwa mem options)
 * -dbq 0: in the (rare) case that some reads have missing qualities for certain bases, default the quality to 0
 * -nofix: GATK will complain and terminate if some reads' mapping quality is over ~60. This option will let GATK keep going.

# Discovering Context-Specific Sequencing Errors

Certain sequence contexts induce errors in next-generation sequencing reads, as detailed in our publications:

Manuel Allhoff, Alexander Schoenhuth, Marcel Martin, Ivan G. Costa, Sven Rahmann and Tobias Marschall. Discovering motifs that induce sequencing errors. BMC Bioinformatics (proceedings of RECOMB-seq), 2013, 14(Suppl 5):S1, DOI: [10.1186/1471-2105-14-S5-S1](http://dx.doi.org/10.1186/1471-2105-14-S5-S1).

On this page, we maintain the source code of our program to discover error-causing motifs.

## Dependencies

 * python 2.7 (tested with Python 2.7.3)
 * HTSeq
 * pysam
 * scipy
 * rpy2 

## Pipeline Dependencies TODO
## Cluster Dependencies

 * python 2.6 (tested with Python 2.6.6)
 * HTSeq
 * pysam
 * scipy 0.15
 * numpy > 1.5.1
 * rpy2 

## Installation

 1. install required python packages (see Dependencies)
 2. checkout the code from the [git repository](https://github.com/zbarni/costalab.git)

## Analysis

Run

```
discovering_cse.py -h
```

to show the help information.

Our tool considers exactly one chromosome in the genome for the analysis. The default chromosome is 'chr1'. You can change it with the option `-c`. If the genome does not have chromosomes (for example E. Coli or B. Subtilis genomes), you do not have to use the `-c` option.

We run the tool, for instance, with the following command:

```
discovering_cse.py hg19.fasta experiment.bam 6 1 -d 0 -c chr10 > results.data
```

Here, we consider chromosome 10 (`-c chr10`) of the human genome (`hg19.fasta`) and search for 6-grams with one allowed N. The analysis is based on the aligned reads which are contained in 'experiment.bam'. Moreover, we do not filter the output (`-d 0`) which is stored in 'results.data'

The output `results.data` is a tab delimited text file and looks like:


| Sequence | Occurrence | Forward Match | Backward Match | Forward Mismatch | Backward Mismatch | Strand Bias Score | FER (Forward Error Rate) | RER (Reverse Error Rate) | ERD (Error rate Difference) |
|----------|------------|---------------|----------------|------------------|-------------------|-------------------|--------------------------|--------------------------|-----------------------------|
| CCANTC   |      12384 |            35 |             37 |             331  |                9  |    12.3345029322  |          0.894557485622  |          0.222222222222  |              0.520215753219 |

We obtain the 6-gram `CCANTC` wich occurs `12384` times in the genome. It gives the following 2x2 contingency table:

|          | Match | Mismatch |
|----------|-------|----------|
|   Foward |    35 |      331 |
| Backward |    37 |        9 |

This table corresponds to a strand bias score of `12.3345029322`.
