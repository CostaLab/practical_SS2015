# costalab
Bioinformatics Lab RWTH Aachen SS '15

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
discovering_cse.py hg19.fasta experiment.bam 6 1 -a 1 -c chr10 > results.data
```

Here, we consider chromosome 10 (`-c chr10`) of the human genome (`hg19.fasta`) and search for 6-grams with one allowed N. The analysis is based on the aligned reads which are contained in 'experiment.bam'. Moreover, we do not filter the output (`-a 1`) which is stored in 'results.data'

The output `results.data` is a tab delimited text file and looks like:


| Sequence | Occurrence | Forward Match | Backward Match | Forward Mismatch | Backward Mismatch | Strand Bias Score | FER (Forward Error Rate) | RER (Reverse Error Rate) | ERD (Error rate Difference) |
|----------|------------|---------------|----------------|------------------|-------------------|-------------------|--------------------------|--------------------------|-----------------------------|
| CCANTC   |      12384 |            35 |             37 |             331  |                9  |    12.3345029322  |          0.894557485622  |          0.222222222222  |                520215753219 |

We obtain the 6-gram `CCANTC` wich occurs `12384` times in the genome. It gives the following 2x2 contingency table:

|          | Match | Mismatch |
|----------|-------|----------|
|   Foward |    35 |      331 |
| Backward |    37 |        9 |

This table corresponds to a strand bias score of `12.3345029322`.
