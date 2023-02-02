#!/bin/bash

workdir=$1
outfile=$2
errfile=$3
script=$4
echo "cd ${workdir} ; bash ${script} > ${outfile} 2> ${errfile}" > /tmp/dbTaskPipe_2024425 
echo "Submitted batch job $RANDOM$RANDOM"
