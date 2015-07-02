#!/bin/bash

if [ ! -d "./organisms" ]
then
	echo "You must cd into the data_analysis subdir, ie you must see the 'organisms' directory"
	exit 1
fi

cd organisms

for f in */*
do
	cd $f

	org_name=`dirname $f`
	org_strain=`basename $f`

	echo
	echo "Summarising stats for $org_name $org_strain"

	BAMS=`find . -name *.bam`

	for bam in BAMS
	do
	done
done