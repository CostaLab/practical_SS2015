#!/bin/bash

########################################################################
# Script implementing a typical bioinformatics pipeline by calling
# external tools on a reference genome and a reads file,
# and producing alignment as well as calling SNPs.
#
# Copyright (C) 2015  Fabio Ticconi
# based on a pipeline from Manuel Allhoff
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

set -o pipefail

# early fail if programs not installed
for c in bwa sra-stat fastq-dump samtools picard-tools gatk fastqc
do
    command -v $c >/dev/null 2>&1 || { echo >&2 "I require \"$c\" (named in this exact way!) but it's not installed. Aborting."; exit 1; }
done

## start option parsing
unset FASTA
unset SRA

# defaults
PLOIDY=1
GATKOPT="-filterNoBases -rf BadCigar"
GATKOPT2=""
MEM=false
SW=false
BWAOPT="-t 4"
FORCE_SINGLE=false
JOIN_FASTQ=false
MINQ="0"
SRANOCHECK=false
LOG="output.log"
SKIP_FASTQC=false

if [ $# == 0 ]
then
    echo `basename $0` "-ref sequence.fasta -sra reads.sra [-p N] [-fix] [-mem] [-mem-pacbio] [-fs] [-dbq N] [-minQ N]"
    exit 1
fi

echo "######################" | tee -a $LOG
echo "## Starting pipeline #" | tee $LOG
echo `basename $0` $@ | tee -a $LOG
echo -e "######################\n" | tee -a $LOG

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -p|--ploidy)
            PLOIDY="$2"
            shift
            ;;
        -ref|--fasta-reference)
            FASTA="$2"
            shift
            ;;
        -sra|--sra-reads)
            SRA="$2"
            shift
            ;;
        -sra-nocheck|--dont-check-sra-file-exists)
            # we need the "SRA" file to extract the READS name,
            # but it doesn't have to exist (assuming we already have the fastq files)
            SRANOCHECK=true
            ;;
        -fix|--fix-qualities)
            GATKOPT=${GATKOPT}" --fix_misencoded_quality_scores"
            ;;
        -nofix|--allow-bad-qualities)
            GATKOPT=${GATKOPT}" --allow_potentially_misencoded_quality_scores"
            ;;
        -mem|--use-bwa-mem)
            MEM=true
            ;;
        -mem-pacbio|--use-bwa-mem-for-pacbio)
            MEM=true
            BWAOPT=${BWAOPT}" -x pacbio"
            ;;
        -mem-ont2d|--use-bwa-mem-for-oxford-nanopore)
            MEM=true
            BWAOPT=${BWAOPT}" -x ont2d"
            ;;
        -sw|--use-bwa-sw)
            SW=true
            ;;
        -fs|--force-single-reads)
            # if there is a single-reads .fastq file, use only that one.
            # otherwise, if there are two, remove the first (ie, is an adaptor) and use
            # the second one as single reads
            FORCE_SINGLE=true
            ;;
        -fs-join|--force-single-join)
            # similar to -fs, but if there are two .fastq files, it "cats" them together.
            # this is useful for certain SRA files that are technically PAIRED-END,
            # but too dodgy to be used in that way: we force them to be treated as single.
            FORCE_SINGLE=true
            JOIN_FASTQ=true
            ;;
        -dbq|--default-base-qualities)
            GATKOPT2="--defaultBaseQualities $2"
            shift
            ;;
        -minQ|--min-mapping-quality)
            MINQ="$2"
            shift
            ;;
	    -skip-fastqc)
	       SKIP_FASTQC=true
	       ;;
        *)
            # unknown option
            echo "## Unknown option: $key" | tee -a $LOG
            echo "## Terminating script" | tee -a $LOG
            exit 1
            ;;
    esac
    shift
done

if [ -z $FASTA ] || [ ! -f $FASTA ]
then
    echo "## Error: must provide a FASTA reference" | tee -a $LOG
    exit 1
elif [ -z $SRA ] || [ ! -f $SRA ]
then
    if [ $SRANOCHECK == false ]
    then
        echo "## Error: must provide a SRA file (use -sra together with -sra-nocheck if you already have .fastq files)" | tee -a $LOG
        exit 1
    fi
fi
## end option parsing

if [ $SRANOCHECK == false ]
then
    echo "## Saving SRA statistics.." | tee -a $LOG
    # generate SRA statistics
    sra-stat --xml -s ${SRA} > `basename ${SRA}`.stats || exit 1
fi

# reference file in FASTA format
REF=`basename ${FASTA} .fasta`

# reads file in SRA format
READS=`basename ${SRA} .sra`

FASTQ1=${READS}_1.fastq
FASTQ2=${READS}_2.fastq
FASTQ3=${READS}_3.fastq
FASTQ4=${READS}_4.fastq
FASTQ=${READS}.fastq

if [ -f $FASTQ ] && [ $FORCE_SINGLE == true ]
then
    echo "## Forcing use of only one fastq file: $FASTQ" | tee -a $LOG

    # fastq files will not be extracted from SRA,
    # and alignment will proceed as from a non-paired
    # dataset

    rm $FASTQ1 &>/dev/null
    rm $FASTQ2 &>/dev/null
    rm $FASTQ3 &>/dev/null
    rm $FASTQ4 &>/dev/null

    FASTQ1=$FASTQ
fi

# We should skip the fastq-dump if $FASTQ1 or $FASTQ2
# exist, but absolutely repeat it if $FASTQ3 or $FASTQ4 exist.
if [ ! -f $FASTQ1 ] || [ -f $FASTQ3 ] || [ -f $FASTQ4 ]
then
    echo "## Dumping fastq files from $SRA" | tee -a $LOG

    # creates fastq split-files
    fastq-dump --split-files $SRA |& tee -a $LOG || exit 1

    if [ -f $FASTQ4 ]
    then
        echo "## Four fastq files generated: removing 1st and 3rd" | tee -a $LOG
        # 4 fastq files means:
        # 1. adaptor
        # 2. real
        # 3. adaptor
        # 4. real
        mv $FASTQ2 $FASTQ1 |& tee -a $LOG || exit 1
        mv $FASTQ4 $FASTQ2 |& tee -a $LOG || exit 1

        rm $FASTQ3 |& tee -a $LOG || exit 1
    elif [ -f $FASTQ3 ]
    then
        echo "## Three fastq files generated: removing 2nd" | tee -a $LOG
        # 3 fastq files means:
        # 1. real
        # 2. adaptor
        # 3. real
        mv $FASTQ3 $FASTQ2 |& tee -a $LOG || exit 1
    elif [ -f $FASTQ2 ]
    then
        echo "## Two fastq files generated" | tee -a $LOG
    else
        echo "## One fastq file generated" | tee -a $LOG
    fi

    if [ $FORCE_SINGLE == true ]
    then
        if [ -f $FASTQ2 ]
        then
            if [ $JOIN_FASTQ == true ]
            then
                echo "## Joining the two paired fastq files into a single file: $FASTQ" | tee -a $LOG
                cat $FASTQ1 $FASTQ2 > $FASTQ || exit 1
                rm $FASTQ1 $FASTQ2 |& tee -a $LOG || exit 1
            else
                echo "## Removing 1st fastq, using only 2nd" | tee -a $LOG
                rm $FASTQ1 |& tee -a $LOG || exit 1
                mv $FASTQ2 $FASTQ |& tee -a $LOG || exit 1
            fi
        else
            echo "## Renaming $FASTQ1 into $FASTQ" | tee -a $LOG
            mv $FASTQ1 $FASTQ |& tee -a $LOG || exit 1
        fi

        FASTQ1=$FASTQ
    fi
fi

echo "## Creating bwa index for ${FASTA}" | tee -a $LOG

bwa index $FASTA |& tee -a $LOG || exit 1

# generate BAM
if [ $MEM == true ]
then
    echo "## Generating sam file with bwa-mem aligner" | tee -a $LOG
    bwa mem $BWAOPT $FASTA *.fastq 2>&1 > tmp.sam | tee -a $LOG || exit 1
elif [ $SW == true ]
then
    echo "## Generating sam file with bwa-sw aligner" | tee -a $LOG
    bwa bwasw $BWAOPT $FASTA *.fastq 2>&1 > tmp.sam | tee -a $LOG || exit 1
else
    echo "## Generating sam file with bwa-backtrack aligner" | tee -a $LOG
    bwa aln $BWAOPT $FASTA $FASTQ1 2>&1 > tmp1.sai | tee -a $LOG || exit 1

    if [ -f $FASTQ2 ]
    then
        bwa aln $BWAOPT $FASTA $FASTQ2 2>&1 > tmp2.sai | tee -a $LOG || exit 1

        bwa sampe $FASTA tmp1.sai tmp2.sai $FASTQ1 $FASTQ2 2>&1 > tmp.sam | tee -a $LOG || exit 1
    else
        bwa samse $FASTA tmp1.sai $FASTQ1 > tmp.sam 2>&1 | tee -a $LOG || exit 1
    fi
fi

echo "## Converting sam file to BAM, sorting and indexing" | tee -a $LOG
samtools view -hbSq $MINQ tmp.sam 2>&1 > tmp.bam | tee -a $LOG || exit 1
samtools sort tmp.bam tmp.sorted |& tee -a $LOG || exit 1
samtools index tmp.sorted.bam |& tee -a $LOG || exit 1

# remove duplicates
echo "## Removing duplicates" | tee -a $LOG
samtools rmdup tmp.sorted.bam tmp.rmdup.sorted.bam |& tee -a $LOG || exit 1

# add readgroup
echo "## Adding read groups and indexing" | tee -a $LOG
picard-tools AddOrReplaceReadGroups INPUT=tmp.rmdup.sorted.bam OUTPUT=tmp.addrg.rmdup.sorted.bam RGLB="rglib" RGPL="rgpl" RGPU="rgpu" RGSM="rgsm" VALIDATION_STRINGENCY=SILENT |& tee -a $LOG || exit 1
samtools index tmp.addrg.rmdup.sorted.bam |& tee -a $LOG || exit 1

# realign near indels
echo "## Realigning near indels and re-indexing" | tee -a $LOG
picard-tools CreateSequenceDictionary R=$FASTA O=${REF}.dict |& tee -a $LOG || exit 1
samtools faidx $FASTA |& tee -a $LOG || exit 1

gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T RealignerTargetCreator -o help.intervals $GATKOPT |& tee -a $LOG || exit 1
gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T IndelRealigner -targetIntervals help.intervals -o ${READS}.bam $GATKOPT $GATKOPT2 |& tee -a $LOG || exit 1
samtools index ${READS}.bam |& tee -a $LOG || exit 1

# generate stats for clean BAM
echo | tee -a $LOG
echo "##################" | tee ${READS}.bam.stats    | tee -a $LOG
echo "# Reads coverage #" | tee -a ${READS}.bam.stats $LOG
echo "##################" | tee -a ${READS}.bam.stats $LOG
if [ $MINQ == 0 ]; then MAPQ=1; else MAPQ=$MINQ; fi
echo "(only counts reads where MAPQ >= ${MAPQ})" | tee -a ${READS}.bam.stats $LOG
GENOME_LENGTH=`samtools view -H ${READS}.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`
echo "Genome length: ${GENOME_LENGTH}" | tee -a ${READS}.bam.stats $LOG
samtools depth -Q${MAPQ} ${READS}.bam | awk -v glen="$GENOME_LENGTH" '{if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; sum+=$3; sumsq+=$3*$3} END { print "min = ",min; print "max = ",max; print "% cov. = ",(NR/glen)*100; print "Avg. base cov. = ",sum/glen; print "Stdev base cov. = ",sqrt(sumsq/glen - (sum/glen)**2)}' | tee -a ${READS}.bam.stats $LOG
echo | tee -a ${READS}.bam.stats $LOG

echo "Reads length (MAPQ >= ${MAPQ}):" | tee -a ${READS}.bam.stats $LOG
echo "min-max avg stdev" | tee -a ${READS}.bam.stats $LOG
samtools view -q $MAPQ ${READS}.bam | awk '{len=length($10); if(max==""){max=min=len}; if(max<len){max=len}; if (min>len){min=len}; sum+=len; sumsq+=len*len} END {avg=sum/NR;sdev=sqrt(sumsq/NR - avg**2);print min"-"max" "avg" +- "sdev}' | tee -a ${READS}.bam.stats $LOG
echo | tee -a ${READS}.bam.stats $LOG

echo "flagstats:"                       | tee -a ${READS}.bam.stats $LOG
samtools flagstat ${READS}.bam         |& tee -a ${READS}.bam.stats $LOG
echo -e "\nidxstats:"                   | tee -a ${READS}.bam.stats $LOG
echo -e "Chr\tlength\tmapped\tunmapped" | tee -a ${READS}.bam.stats $LOG
samtools idxstats ${READS}.bam         |& tee -a ${READS}.bam.stats $LOG
echo | tee -a ${READS}.bam.stats $LOG

if [ $SKIP_FASTQC == false ]
then
    # generate other statistics
    echo "## Running fastqc for web-based statistics" | tee -a $LOG
    fastqc ${READS}.bam |& tee -a $LOG
fi

# SNP calling
echo "## SNP calling with GATK" | tee -a $LOG
gatk -ploidy $PLOIDY -I ${READS}.bam -R $FASTA -T UnifiedGenotyper -o ${READS}-snps.vcf $GATKOPT |& tee -a $LOG || exit 1

echo "SNPs found: "`egrep -c "^[^#]" *.vcf` |& tee -a $LOG
