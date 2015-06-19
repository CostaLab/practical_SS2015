#!/usr/bin/env bash

if [[ $# < 1 ]]
then
    echo "`basename $0` [cat/bsub] command args"
    exit -1
fi

command -v $2 >/dev/null 2>&1 || { echo >&2 "Command '$1' not available"; exit -1; }

bin=`which $2`
args="${*:3}"

$1 <<EOF
#!/usr/bin/env zsh
### Job name
#BSUB -J $2
 
### File / path where STDOUT & STDERR will be written
### %J is the job ID, %I is the array ID
#BSUB -o $2.log.%J.%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 3:0
 
### Request memory you need for your job in TOTAL in MB
#BSUB -M 15000
 
### Change to the work directory
C=\$LSB_JOBINDEX
cd $PWD

source /home/${USER}/.zshrc

### Execute your application
$bin $args
EOF

