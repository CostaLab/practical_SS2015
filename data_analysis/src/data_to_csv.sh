#!/bin/bash

mkdir -p csv

for i in `ls results_*`
do
	echo "Motif,Occurrences,StrandBiasScore,ERD" > csv/`basename $i .data`.csv
	cat $i | egrep "^[^ #]" | awk '{print $1 "," $2 "," $7 "," $10}' >> csv/`basename $i .data`.csv
done
