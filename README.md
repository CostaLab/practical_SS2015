# Bioinformatics Lab RWTH Aachen SS '15

# Align, statistics, SNP calling

The script bio_pipeline.sh allows an easy and almost-painless alignment procedure, as well as statistics generation and finally SNP calling via GATK.

At the very least, it requires two files:

 * a FASTA reference genome (option -ref) (this MUST be in the current directory)
 * a SRA reads file (option -sra) (this can be in any directory, or even "fake" with option -sra-nocheck)

```
bio_pipeline.sh -ref sequence.fasta -sra reads.sra
```

Without other options, the script will:

 * assume the organism is haploid (ie, only has one set of chromosomes)
 * extract the SRA file to FASTQ files, and remove (in simple cases) the adaptor-only files
  * the previous step also infers whether the SRA is paired-ends or single-ends
 * align the reads to the reference genome, using bwa-backtrack with default options
 * convert SAM file to BAM
 * remove duplicates, re-align near indels, produce a clean BAM file
 * generate plenty of statistics for the final BAM file
 * run GATK SNP calling and print the final number of SNPs

Note: the name of the SRA file determines the prefix for the fastq files, the vcf file and so on.

All the files are kept: it is up to you to delete the temporary files, or the unnecessary ones. A log is generated to output.log, containing both the stdout and stderr produced during the execution of bio_pipeline.sh.

The bare command presented above is OK for Illumina short-reads, either paired or single ends.

Among the options, some are notable:

 * -mem: use bwa-mem instead of backtrack (very important for Ion Torrent, 454)
 * -mem-pacbio: a special option for PacBio, instead of the previous one (adds "-x pacbio" to the bwa mem options)
 * -mem-ont2d: a special option for Oxford Nanopore reads (it simply adds "-x ont2d" to the bwa mem options)
 * -dbq 0: in the (rare) case that some reads have missing qualities for certain bases, default the quality of the missing bases to 0 instead of ditching the read. Might be useful for 454 reads
 * -nofix: GATK will complain and terminate if some reads' mapping quality is over ~60. This option will let GATK keep going. It is advised that you first check if GATK was right to complain (eg, the quality encoding is non-standard)

## Examples

### Remove technical reads

Many SRA files contain multiple sequences, all dumped when using ```fastq-dump --split-files reads.sra```. The cases we have experienced are:

1. 4 files: it means reads_1.fastq is technical, reads_2.fastq is not technical, reads_3.fastq is technical, reads_4.fastq is not technical. The first and third file must be removed, the second and fourth must be used in a paired-end alignment
2. 3 files: it means reads_1.fastq is not technical, reads_2.fastq is technical, reads_3.fastq is not technical. The second file must be removed and the other two be used in a paired end alignment
3. 2 files: either they are both non technical, or reads_1.fastq is technical while reads_2.fastq is not. In the first case, both must be used in paired end alignment. In the second case, reads_1.fastq must be removed and reads_2.fastq be used in a single-end alignment

Case 1 and 2 are automatically handled. Case 3, by default, assume the two reads are biological reads and use them in a paired-end alignment. However, if we use the option -fs, it will remove reads_1.fastq and rename reads_2.fastq into reads.fastq, then using it in a single-end alignment.

```
bio_pipeline.sh -ref sequence.fasta -sra reads.sra -fs
# the only fastq file present here is reads.fastq
```

### FASTQ files already present

Suppose you are directly using fastq files instead of SRA files. What you want is to skip the SRA file check or the script won't run. *You need to always have either one fastq or two fastq files, never less or more*.

```
bio_pipeline.sh -ref sequence.fasta -sra reads.sra -sra-nocheck
```

Note how you still need to specify the sra file name. This is used to extract the prefix (in this case, "reads") so as to fetch the correct fastq files.

```
# The above command will look, at least, for:
reads_1.fastq
# If the following is also present, it will proceed with
# a paired-end alignment. Otherwise, it will only use the first
# and proceed to single-end alignment
reads_2.fastq
```

If you have single-end fastq files in the more compact format "reads.fastq", you just need to add the option -fs.

```
bio_pipeline.sh -ref sequence.fasta -sra reads.sra -sra-nocheck -fs
```

This will look for reads.fastq and, if found, proceed with a single-end alignment.

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
