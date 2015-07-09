#!/bin/bash

# Goes inside every organism/platform combination,
# then it takes all available .data files (in all three possible levels:
# consensus, experiment and run) and merges them together.
# this means that all the contingency tables will be merged and nothing
# will be lost

if [ ! -d "./organisms" ]
then
	echo "You must cd into the data_analysis subdir, ie you must see the 'organisms' directory"
	exit 1
fi

restype="results"

cd organisms

base="$PWD"
for f in */*
do
	cd $f

	org_name=`dirname $f`
	org_strain=`basename $f`

	echo
	echo "# Merging motifs for $org_name $org_strain #"
	echo

	org_base="$PWD"
	shopt -s nullglob
	for platform in `ls *{GAII,HiSeq,PGM_mem,MiSeq,PacBio,MinION}* -d`
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
		elif [[ $platform == MinION ]]
		then
			platform=minion
		fi

		# type=commonstrict disabled for now
		for type in merged
		do
			if [ $type == "commonstrict" ]
			then
				opts=""
			else
				opts="-r simple"
			fi
			
			for d in 0 0.05
			do
				d2=`echo $d | tr -d .`
				
				n=2
				for q in 4 8 10
				do
					files=`shopt -s nullglob; echo {$restype*,*/$restype,*/*/$restype}/results_${q}-grams_${n}n.data`
					simple_motif_merger.py $q $n $files -d $d $opts > ${platform}_${type}_${restype}_${q}-grams_${n}n_d${d2}.data
				done
				files=`shopt -s nullglob; echo {$restype*,*/$restype,*/*/$restype}/results_8-grams_4n.data`
				simple_motif_merger.py 8 4 $files -d $d $opts > ${platform}_${type}_${restype}_8-grams_4n_d${d2}.data
			done
		done

		cd $org_base
	done
	shopt -u nullglob

	cd $base
done

