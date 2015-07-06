#!/bin/bash

# Indel calling

base="$PWD"
for f in `find . -name "*.bam" | grep -v -E "(staph|454_|_old|_sw)"`
do
    echo $f
    
    dir=`dirname $f`
    cd $dir

    bam=`basename $f`

    name=`basename $f .bam`

    if [ -f ${name}-indels.vcf ]
    then
        echo "indels vcf already present, skipping"
        cd $base
        continue
    fi

    if [ ! -f ${bam}.bai ]
    then
        echo "generating bam index.."
        samtools index $bam
    fi

    name=`basename $f .bam`

    if [ -f ../*.fasta ]
    then
        cp ../*.fasta .
    elif [ -f ../../*.fasta ]
    then
        cp ../../*.fasta .
    elif [ -f ../../../*.fasta ]
    then
        cp ../../../*.fasta .
    else
        echo $dir
        echo $bam
        echo $name
        echo "missing fasta file for $f!"
        exit 1
    fi

    FASTA=`basename *.fasta .fasta`.fasta
    samtools faidx $FASTA
    picard-tools CreateSequenceDictionary R=$FASTA O=`basename $FASTA .fasta`.dict
    
    gatk -ploidy 1 -I $bam -R $FASTA -T UnifiedGenotyper -o ${name}-indels.vcf -glm INDEL -rf BadCigar -filterNoBases --allow_potentially_misencoded_quality_scores --defaultBaseQualities 0

    rm *.fasta
    rm *.dict
    rm *.fai
    
    cd $base
done
