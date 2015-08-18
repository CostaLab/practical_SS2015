#!/usr/bin/env bash

if (( $# < 6 ))
then
    echo "`basename $0` [cat/bsub] proj name W M command args"
    exit -1
fi

command -v $6 >/dev/null 2>&1 || { echo >&2 "Command '$6' not available"; exit -1; }

bin=`which $6`
args="${*:7}"

proj="$2"
name="$3"
W="$4"
M="$5"

$1 <<EOF
#!/usr/bin/env zsh
### Job name
#BSUB -J $name

##BSUB -P $proj
 
### File / path where STDOUT & STDERR will be written
### %J is the job ID, %I is the array ID
#BSUB -o ${name}.log.%J.%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W $W
 
### Request memory you need for your job in TOTAL in MB
#BSUB -M $M
 
### Change to the work directory
C=\$LSB_JOBINDEX
cd $PWD

source /home/${USER}/.zshrc

### Execute your application
$bin $args
EOF
