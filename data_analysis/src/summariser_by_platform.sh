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

TMP="${PWD}/.tmp_`basename $0 .sh`"
rm $TMP.* &> /dev/null

echo \
"Platform,Organisms,BAMs,WAvg Reads,PStdDev Reads,WAvg Read Depth,PStdDev Read Depth,WAvg Read Length,\
PStdDev Read Length,Avg SNPs,StdDev SNPs,Avg Indels,StdDev Indels,Tot SNPs-motifs,Tot SNPs-motifs 0.05" | tee $FILE

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
	# bnum: number of bam files considered

	org_base="$PWD"
	shopt -s nullglob
	for platform in `ls *{GAII,HiSeq,PGM_mem,MiSeq,PacBio,MinIon}* -d`
	do
		rm $TMP.* &> /dev/null
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

		total_motifs_file="${platform}_merged_results_${q}-grams_${n}n_d0.data"
		tot_motifs=`egrep -c "^[^#]" $total_motifs_file`

		total_motifs005_file="${platform}_merged_results_${q}-grams_${n}n_d005.data"
		tot_motifs005=`egrep -c "^[^#]" $total_motifs005_file`

		BAMS=`find . -name "*.bam"`
		bnum=`ls $BAMS | wc -l`

		for bam in $BAMS
		do
			echo $bam

			# we accumulate the bases covered, so that at the end we will
			# have the total breadth of coverage for all bams
			# (it will have to be uniquely sorted to remove duplicates)
			samtools depth -Q0 $bam > ${TMP}.hcov.new

			if [[ `cat ${TMP}.hcov.new | wc -l` == 0 ]]
			then
				# if the bam doesn't have mapped reads,
				# we don't include it
				bnum=$((bnum - 1))
				continue
			fi

			cat ${TMP}.hcov.new | awk '{print $1" "$2}' >> ${TMP}.hcov

			# we then re-use the extracted read depth data from two commands above
			res=`cat ${TMP}.hcov.new | awk -v glen="$glen" '{sum+=$3;sumsq+=$3*$3} END {print "rdsize="NR;print "rdavg="sum/glen;print "rdstd="sqrt(sumsq/glen - (sum/glen)**2)}'`
			for r in $res
			do
				eval $r
			done

			res2=`samtools view -F 4 $bam | awk '{len=length($10);sum+=len;sumsq+=len*len} END {avg=sum/NR;sdev=sqrt(sumsq/NR - avg**2);print "rlavg="avg,"rlstd="sdev;print "rtot="NR}'`
			for r in $res2
			do
				eval $r
			done

			# rdsize: number of bases covered by at least 1 read
			# rdavg, rdstd: read depth
			# rtot: number of reads
			# rlavg, rlstd: read length

			dir=`dirname $bam`
			name=`basename $bam .bam`
			snps_count=`egrep -c "^[^#]" ${dir}/${name}-snps.vcf`
			indel_count=`egrep -c "^[^#]" ${dir}/${name}-indels.vcf`

            # for read depth and read length, we need to calculate the combined variance.
            # we need to accumulate sample size, avg and std
            echo "$rdsize $rdavg $rdstd" >> ${TMP}.rd
            echo "$rtot $rlavg $rlstd" >> ${TMP}.rl
            echo "$snps_count $indel_count" >> ${TMP}.snpind

		done

		# number of bases covered by at least 1 read, in at least 1 bam,
		# over the length of the genome. Should be 100% most of the times
		hcov=`cat ${TMP}.hcov | sort -k 1,2 -u | wc -l`
		hcovperc=`echo "scale=8; ($hcov / $glen)*100" | bc -l`

		res=`cat ${TMP}.rd | awk '{wavgi+=($1*$2);sizes+=$1;pvar+=(($1 - 1)*$3*$3)} END {print "wavgrd="(wavgi/sizes);print "pstdrd="sqrt(pvar/(sizes - NR))}'`
		for r in $res
		do
			eval $r
		done
		# wavgrd: weighted average of read depths
		# pstdrd: pooled standard deviation of read depths

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

		echo $org_name,$org_strain,$glen,$gc_stats,$max_kmers,$cov_kmers,$per_kmers,$platform,$bnum,$hcovperc,$rtotavg,$rtotstd,$wavgrd,$pstdrd,$wavgrl,$pstdrl,$snpavg,$snpstd,$indavg,$indstd,$tot_motifs,$tot_motifs005 | tee -a $FILE

		cd $org_base
	done
	shopt -u nullglob

	cd $base
done
