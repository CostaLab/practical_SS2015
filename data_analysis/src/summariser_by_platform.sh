#!/bin/bash

if [ ! -d "./organisms" ]
then
	echo "You must cd into the data_analysis subdir, ie you must see the 'organisms' directory"
	exit 1
fi

if [[ $# < 2 ]]
then
    echo "`basename $0` Q N"
    exit 1
fi

q=$1
n=$2

FILE="${PWD}/stats_by_platform_${q}_${n}.csv"

DIR="${PWD}/final_results/bs_bp_ec_sa_se_vc_pf"

TMP="${PWD}/.tmp_`basename $0 .sh`"
rm $TMP.* &> /dev/null

echo \
"Platform,Organisms,BAMs,WAvg Reads,PStdDev Reads,WAvg Read Length,\
PStdDev Read Length,Avg SNPs,StdDev SNPs,Avg Indels,StdDev Indels,\
Tot SNPs-motifs,Tot SNPs-motifs 0.05" | tee $FILE

cd organisms

for platform in GAII HiSeq PGM_mem MiSeq PacBio MinIon
do
	if [[ $platform == GAII ]]
	then
		platform2=gaii-x
	elif [[ $platform == HiSeq ]]
	then
		platform2=hiseq
	elif [[ $platform == MiSeq ]]
	then
		platform2=miseq
	elif [[ $platform == PGM_mem ]]
	then
		platform2=ion-torrent
	elif [[ $platform == PacBio ]]
	then
		platform2=pacbio
	elif [[ $platform == MinIon ]]
	then
		platform2=minion
	fi

	echo
	echo "# Summarising stats for $platform2 #"
	echo

	orgnum=0
	bnumtot=0
	rm $TMP.* &> /dev/null

	shopt -s nullglob
	base="$PWD"
	for f in `ls -d */*/*${platform}*`
	do
		cd $f

		orgnum=$((orgnum + 1))

		total_motifs_file="${DIR}/merged_from_merged/${platform2}_merged_results_${q}-grams_${n}n_d"
		tot_motifs=`egrep -c "^[^#]" ${total_motifs_file}0.data`
		tot_motifs005=`egrep -c "^[^#]" ${total_motifs_file}005.data`

		common_motifs_file="${DIR}/commonstrict_from_merged/${platform2}_commonstrict_results_${q}-grams_${n}n_d"
		common_motifs=`egrep -c "^[^#]" ${common_motifs_file}0.data`
		common_motifs005=`egrep -c "^[^#]" ${common_motifs_file}005.data`

		BAMS=`find . -name "*.bam"`
		bnum=`ls $BAMS | wc -l`

		for bam in $BAMS
		do
			echo $bam

			res2=`samtools view -F 4 $bam | awk '{len=length($10);sum+=len;sumsq+=len*len} END {avg=sum/NR;sdev=sqrt(sumsq/NR - avg**2);print "rlavg="avg,"rlstd="sdev;print "rtot="NR}'`
			for r in $res2
			do
				eval $r
			done

			# rtot: number of reads
			# rlavg, rlstd: read length

			dir=`dirname $bam`
			name=`basename $bam .bam`
			snps_count=`egrep -c "^[^#]" ${dir}/${name}-snps.vcf`
			indel_count=`egrep -c "^[^#]" ${dir}/${name}-indels.vcf`

            # for read length, we need to calculate the combined variance.
            # we need to accumulate sample size, avg and std
            echo "$rtot $rlavg $rlstd" >> ${TMP}.rl
            echo "$snps_count $indel_count" >> ${TMP}.snpind
		done

		cd $base
	done
	shopt -u nullglob

	bnumtot=$((bnumtot + bnum))

	res=`cat ${TMP}.rl | awk '{sum+=$1;sumsq+=$1*$1;wavgi+=$1*$2;sizes+=$1;pvar+=(($1 - 1)*$3*$3)} END {avg=sum/NR;print "rtotavg="avg; print "rtotstd="sqrt(sumsq/NR - avg**2);print "wavgrl="(wavgi/sizes);print "pstdrl="sqrt(pvar/(sizes - NR))}'`
	for r in $res
	do
		eval $r
	done
	# rtotavg: average of number of reads
	# rtotstd: standard deviation of number of reads
	# wavgrl: weighted average of read lengths
	# pstdrl: pooled standard deviation of read lengths

	res=`cat ${TMP}.snpind | awk '{snpsum+=$1;snpsumsq+=$1*$1;indsum+=$2;indsumsq+=$2*$2} END {snpavg=snpsum/NR;indavg=indsum/NR;print "snpavg="snpavg;print "snpstd="sqrt(snpsumsq/NR - snpavg**2);print "indavg="indavg;print "indstd="sqrt(indsumsq/NR - indavg**2)}'`
	for r in $res
	do
		eval $r
	done
	# snpavg: average of snp count
	# snpstd: st. dev. of snp count
	# indavg: average of indel count
	# indstd: st. dev. of indel count

	echo $platform2,$orgnum,$bnumtot,$rtotavg,$rtotstd,$wavgrl,$pstdrl,$snpavg,$snpstd,$indavg,$indstd,$tot_motifs,$tot_motifs005,$common_motifs,$common_motifs005 | tee -a $FILE
done
