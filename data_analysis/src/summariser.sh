#!/bin/bash

if [ ! -d "./organisms" ]
then
	echo "You must cd into the data_analysis subdir, ie you must see the 'organisms' directory"
	exit 1
fi

FILE="${PWD}/stats.csv"

echo "Organism,Strain,Genome Length,%GC,%CpG,CpG Obs/Exp,Q,N,Max Kmers,Covered Kmers,%Cov Kmers,Platform,Tot Motifs,Common Motifs,JC,BAM,Tot Reads,HCov,Avg Read Depth,Avg Read Length,SNPs,Indels" | tee $FILE

q=10
n=2

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

	#bnum=`echo $BAMS | wc -l`

	# glen: genome length
	# bnum: number of bam files considered

	org_base="$PWD"
	shopt -s nullglob
	for platform in `ls *{GAII,HiSeq,PGM_mem,MiSeq,PacBio,MinIon}* -d`
	do
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

		total_motifs_file="${platform}_merged_results_${q}-grams_${n}n_d005.data"
		tot_motifs=`egrep -c "^[^#]" $total_motifs_file`

		common_motifs_file="${platform}_commonstrict_results_${q}-grams_${n}n_d005.data"
		intersection=`egrep -c "^[^#]" $common_motifs_file`

		if [[ $tot_motifs == 0 ]]
		then
			JC=0
		else
			JC=`echo "scale=8; $intersection / $tot_motifs" | bc -l`
		fi

		BAMS=`find . -name "*.bam" | grep -v -E "(staph|454_|_old|_sw)"`

		for bam in $BAMS
		do
			echo $bam
			res=`samtools depth -Q0 $bam | awk -v glen="$glen" '{if(min==""){min=max=$3};if($3>max){max=$3};if($3< min){min=$3};sum+=$3;sumsq+=$3*$3} END {print "hcov="(NR/glen)*100;print "rdmin="min;print "rdmax="max;print "rdavg="sum/glen;print "rdstd="sqrt(sumsq/glen - (sum/glen)**2)}'`
			for r in $res
			do
				eval $r
			done

			res2=`samtools view -F 4 $bam | awk '{len=length($10);if(max==""){max=min=len};if(max<len){max=len};if(min>len){min=len};sum+=len;sumsq+=len*len;rtot+=1} END {avg=sum/NR;sdev=sqrt(sumsq/NR - avg**2);print "rlmin="min,"rlmax="max,"rlavg="avg,"rlstd="sdev;print "rtot="rtot}'`
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

			echo $org_name,$org_strain,$glen,$gc_stats,$q,$n,$max_kmers,$cov_kmers,$per_kmers,$platform,$tot_motifs,$intersection,$JC,$name,$rtot,$hcov,$rdavg,$rlavg,$snps_count,$indel_count | tee -a $FILE
		done

		cd $org_base
	done
	shopt -u nullglob

	cd $base
done