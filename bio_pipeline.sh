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
GATKOPT=""
GATKOPT2=""
MEM=false
BWAOPT="-t 4"
FORCE_SINGLE=false
MINQ=0

if [ $# == 0 ]
then
    echo `basename $0` "-ref sequence.fasta -sra reads.sra [-p N] [-fix] [-mem] [-mem-pacbio] [-fs] [-dbq N] [-minQ N]"
    exit 1
fi

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
            GATKOPT="--fix_misencoded_quality_scores"
            ;;
        -nofix|--allow-bad-qualities)
            GATKOPT="--allow_potentially_misencoded_quality_scores"
            ;;
        -mem|--use-bwa-mem)
            MEM=true
            ;;
        -mem-pacbio|--use-bwa-mem-for-pacbio)
            MEM=true
            BWAOPT=${BWAOPT}" -x pacbio"
            ;;
        -fs|--force-single-reads)
            FORCE_SINGLE=true
            ;;
        -dbq|--default-base-qualities)
            GATKOPT2="--defaultBaseQualities $2"
            shift
            ;;
        -minQ|--min-mapping-quality)
            MINQ="$2"
            shift
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

# generate SRA statistics
sra-stat --xml -s ${SRA} > ${SRA}.stats

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
    echo "# Using only one fastq file: $FASTQ"

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
# exist, but absolutely repeat it if $FASTQ3 or $FASTQ3 exist.
if [ ! -f $FASTQ1 ] || [ -f $FASTQ3 ] || [ -f $FASTQ4 ]
then
    # creates fastq split-files
    fastq-dump --split-files $SRA || exit 1

    if [ -f $FASTQ4 ]
    then
        echo "# Four fastq files generated: removing 1st and 3rd"
        # 4 fastq files means:
        # 1. adaptor
        # 2. real
        # 3. adaptor
        # 4. real
        mv $FASTQ2 $FASTQ1 || exit 1
        mv $FASTQ4 $FASTQ2 || exit 1

        rm $FASTQ3 || exit 1
    elif [ -f $FASTQ3 ]
    then
        echo "# Three fastq files generated: removing 2nd"
        # 3 fastq files means:
        # 1. real
        # 2. adaptor
        # 3. real
        mv $FASTQ3 $FASTQ2 || exit 1

        rm $FASTQ4 &>/dev/null
    elif [ -f $FASTQ2 ]
    then
        echo "# Two fastq files generated"
    else
        echo "# One fastq file generated"
    fi

    if [ $FORCE_SINGLE == true ]
    then
        if [ -f $FASTQ2 ]
        then
            #echo "# Joining the two paired fastq files into a single file: $FASTQ"
            #cat $FASTQ1 $FASTQ2 > $FASTQ || exit 1
            #rm $FASTQ1 $FASTQ2 || exit 1
            echo "# Removing 1st fastq, using only 2nd"
            rm $FASTQ1 || exit 1
            mv $FASTQ2 $FASTQ || exit 1
        else
            echo "# Renaming $FASTQ1 into $FASTQ"
            mv $FASTQ1 $FASTQ || exit 1
        fi

        FASTQ1=$FASTQ
    fi
fi

# generate BAM
if [ $MEM == true ]
then
    bwa index $FASTA || exit 1
    bwa mem $BWAOPT $FASTA *.fastq > tmp.sam || exit 1
else
    bwa index $FASTA || exit 1
    bwa aln $BWAOPT $FASTA $FASTQ1 > tmp1.sai || exit 1

    if [ -f $FASTQ2 ]
    then
        bwa aln $BWAOPT $FASTA $FASTQ2 > tmp2.sai || exit 1

        bwa sampe $FASTA tmp1.sai tmp2.sai $FASTQ1 $FASTQ2 > tmp.sam || exit 1
    else
        bwa samse $FASTA tmp1.sai $FASTQ1 > tmp.sam || exit 1
    fi
fi

samtools view -bSq $MINQ tmp.sam > tmp.bam || exit 1
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

gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T RealignerTargetCreator -o help.intervals $GATKOPT || exit 1
gatk -I tmp.addrg.rmdup.sorted.bam -R $FASTA -T IndelRealigner -targetIntervals help.intervals -o ${READS}.bam $GATKOPT $GATKOPT2 || exit 1
samtools index ${READS}.bam || exit 1

# generate stats for clean BAM
echo "##################" | tee ${READS}.bam.stats
echo "# Reads coverage #" | tee -a ${READS}.bam.stats
echo "##################" | tee -a ${READS}.bam.stats
echo "(excluding MAPQ=0)" | tee -a ${READS}.bam.stats
GENOME_LENGTH=`samtools view -H ${READS}.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`
echo "Genome length: ${GENOME_LENGTH}" | tee -a ${READS}.bam.stats
samtools depth -Q1 ${READS}.bam | awk -v glen="$GENOME_LENGTH" '{if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; sum+=$3; sumsq+=$3*$3} END { print "min = ",min; print "max = ",max; print "% cov. = ",(NR/glen)*100; print "Avg. base cov. = ",sum/glen; print "Stdev base cov. = ",sqrt(sumsq/glen - (sum/glen)**2)}' | tee -a ${READS}.bam.stats
echo | tee -a ${READS}.bam.stats

echo "flagstats:"                         | tee -a ${READS}.bam.stats
samtools flagstat ${READS}.bam            | tee -a ${READS}.bam.stats
echo -e "\nidxstats:"                     | tee -a ${READS}.bam.stats
echo -e "Chr\tlength\tmapped\tunmapped"   | tee -a ${READS}.bam.stats
samtools idxstats ${READS}.bam            | tee -a ${READS}.bam.stats

# generate other statistics
fastqc ${READS}.bam

# SNP calling
gatk -ploidy $PLOIDY -I ${READS}.bam -R $FASTA -T UnifiedGenotyper -o ${READS}-snps.vcf $GATKOPT || exit 1

echo "SNPs found: "`egrep -c "^[^#]" *.vcf`