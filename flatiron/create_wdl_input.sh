#!/bin/bash

c=$1
w=`pwd`

for i in {1..7}
do
	./create_bam_tsv.pl iWES.bam.list ${w}/chr${c}_${i}_tsvs > iWES_chr${c}_${i}_tsv_bams.tsv
done

