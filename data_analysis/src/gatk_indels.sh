#!/bin/bash

# Indel calling

base="$PWD"
for f in `find . -name "*.bam" | grep -v -E "(staph|454|_old|_sw)"`
do
    dir=`dirname $f`
    cd $dir

    bam=`basename $f`
    samtools index $bam

    name=`basename $f .bam`

    if [ -f ../*.fasta ]
    then
        cp ../*.fasta .
    elif [ -f ../../*.fasta ]
    then
        cp ../../*.fasta
    elif [ -f ../../../*.fasta]
    then
        cp ../../../*.fasta
    else
        echo $dir
        echo $bam
        echo $name
        echo "missing fasta file for $f!"
        exit 1
    fi
    
    gatk -ploidy 1 -I $bam -R *.fasta -T UnifiedGenotyper -o ${name}-indels.vcf -glm INDEL -rf BadCigar -filterNoBases --allow_potentially_misencoded_quality_scores --defaultBaseQualities 0
    
    cd $base
done
