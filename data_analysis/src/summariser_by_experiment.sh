#!/bin/bash

if [ ! -d "./organisms" ]
then
	echo "You must cd into the data_analysis subdir, ie you must see the 'organisms' directory"
	exit 1
fi

if [[ $# < 2 ]]
then
    echo "$0 Q N"
    exit 1
fi

q=$1
n=$2

FILE="${PWD}/stats_by_exp_${q}_${n}.csv"

echo "Organism,Strain,Genome Length,%GC,%CpG,CpG Obs/Exp,Max Kmers,Covered Kmers,%Cov Kmers,Platform,BAM,Tot Reads,HCov,Avg Read Depth,StdDev Read Depth,Avg Read Length,StdDev Read Length,SNPs,Indels,SNPs-motifs,SNPs-motifs 0.05" | tee $FILE

cd organisms

base="$PWD"
for f in */*
do
	cd $f

	org_name=`dirname $f`
	org_strain=`basename $f`

	echo
	echo "# Summarising stats for $org_name $org_strain #"
	echo

	glen=`grep "DNA bps" README.md -A 2 | tail -n 1 | awk -F\| '{print $4}' | awk '{print $1}'`

	max_kmers=`grep "## Total ##" ${q}kmers.txt -A 1 | tail -n 1 | awk '{print $7}'`
	cov_kmers=`grep "## Total ##" ${q}kmers.txt -A 1 | tail -n 1 | awk '{print $4}'`
	per_kmers=`grep "## Total ##" ${q}kmers.txt -A 1 | tail -n 1 | awk '{print $8}' | tr -d "()%"`

	gc_stats=`tail -n +2 gc_out.txt | awk '{gc+=$2;cpg+=$3;oer+=$4} END {print (gc/NR)","(cpg/NR)","(oer/NR)}'`

	# glen: genome length

	org_base="$PWD"
	shopt -s nullglob
	for platform in `ls *{GAII,HiSeq,PGM_mem,MiSeq,PacBio,MinIon}* -d`
	do
        echo "Platform: $platform"

		cd $platform

		if [[ $platform == GAII* ]]
		then
			platform=gaii-x
		elif [[ $platform == HiSeq* ]]
		then
			platform=hiseq
		elif [[ $platform == MiSeq* ]]
		then
			platform=miseq
		elif [[ $platform == IonTorr* ]]
		then
			platform=ion-torrent
		elif [[ $platform == PacB* ]]
		then
			platform=pacbio
		elif [[ $platform == MinIon ]]
		then
			platform=minion
		fi

		BAMS=`find . -name "*.bam"`

		for bam in $BAMS
		do
			echo $bam
			res=`samtools depth -Q0 $bam | awk -v glen="$glen" '{sum+=$3;sumsq+=$3*$3} END {print "hcov="(NR/glen)*100,"rdavg="sum/glen,"rdstd="sqrt(sumsq/glen - (sum/glen)**2)}'`
			for r in $res
			do
				eval $r
			done

			res2=`samtools view -F 4 $bam | awk '{len=length($10);sum+=len;sumsq+=len*len;rtot+=1} END {avg=sum/NR;sdev=sqrt(sumsq/NR - avg**2);print "rlavg="avg,"rlstd="sdev,"rtot="rtot}'`
			for r in $res2
			do
				eval $r
			done

			# hcov: horizontal coverage
			# rtot: number of reads
			# rdmin, rdmax, rdavg, rdstd: basic stats of read depth
			# rlmin, rlmax, rlavg, rlstd: basic stats of read length

			dir=`dirname $bam`
			name=`basename $bam .bam`
			snps_count=`egrep -c "^[^#]" ${dir}/${name}-snps.vcf`
			indel_count=`egrep -c "^[^#]" ${dir}/${name}-indels.vcf`
            snps_motifs_count=`egrep -c "^[^ #]" ${dir}/results/results_${q}-grams_${n}n.data`
            snps_motifs_count_005=`egrep "^[^ #]" ${dir}/results/results_${q}-grams_${n}n.data | awk '{if ($10 > 0.05) print $0}' | wc -l`

			echo $org_name,$org_strain,$glen,$gc_stats,$max_kmers,$cov_kmers,$per_kmers,$platform,$name,$rtot,$hcov,$rdavg,$rdstd,$rlavg,$rlstd,$snps_count,$indel_count,$snps_motifs_count,$snps_motifs_count_005 | tee -a $FILE
		done

		cd $org_base
	done
	shopt -u nullglob

	cd $base
done
