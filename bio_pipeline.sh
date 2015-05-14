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

## start option parsing
unset FASTA
unset SRA

# defaults
PLOIDY=1
FIX=""
MEM=false

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
        -fix|--fix-qualities)
            FIX="--fix_misencoded_quality_scores"
            ;;
        -mem|--use-bwa-mem)
            MEM=true
            ;;
        *)
            # unknown option
            ;;
    esac
    shift
done

if [ -z $FASTA ] || [ ! -f $FASTA ]
then
    echo "Must provide a FASTA reference"
    exit
elif [ -z $SRA ] || [ ! -f $SRA ]
then
    echo "Must provide a SRA file"
    exit
fi
## end option parsing

# reference file in FASTA format
REF=`basename ${FASTA} .fasta`

# reads file in SRA format
READS=`basename ${SRA} .sra`

FASTQ1=${READS}_1.fastq
FASTQ2=${READS}_2.fastq

# creates fastq files
fastq-dump --split-files $SRA || exit 1

# generate BAM
if [ $MEM == true ]
then
    bwa index $FASTA || exit 1
    bwa mem -t 4 $FASTA *.fastq > tmp.sam || exit 1
else
    bwa index $FASTA || exit 1
    bwa aln $FASTA $FASTQ1 > tmp1.sai || exit 1

    if [ -f $FASTQ2 ]
    then
        bwa aln $FASTA $FASTQ2 > tmp2.sai || exit 1

        bwa sampe $FASTA tmp1.sai tmp2.sai $FASTQ1 $FASTQ2 > tmp.sam || exit 1
    else
        bwa samse $FASTA tmp1.sai $FASTQ1 > tmp.sam || exit 1
    fi
fi

samtools view -bS tmp.sam > tmp.bam || exit 1
samtools sort tmp.bam tmp.sorted || exit 1
samtools index tmp.sorted.bam || exit 1

# remove duplicates
samtools rmdup tmp.sorted.bam tmp.rmdup.sorted.bam || exit 1

# add readgroup
picard-tools AddOrReplaceReadGroups INPUT=tmp.rmdup.sorted.bam OUTPUT=tmp.addrg.rmdup.sorted.bam RGLB="rglib" RGPL="rgpl" RGPU="rgpu" RGSM="rgsm" VALIDATION_STRINGENCY=SILENT || exit 1
samtools index tmp.addrg.rmdup.sorted.bam || exit 1

# realign near indels
picard-tools CreateSequenceDictionary R=$FASTA O=${REF}.dict || exit 1
samtools faidx $FASTA || exit 1

gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T RealignerTargetCreator -o help.intervals $FIX || exit 1
gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T IndelRealigner -targetIntervals help.intervals -o ${READS}.bam $FIX || exit 1
samtools index ${READS}.bam || exit 1

# SNP calling
gatk -ploidy $PLOIDY -I ${READS}.bam -R $FASTA -T UnifiedGenotyper -o ${READS}-snps.vcf || exit 1

