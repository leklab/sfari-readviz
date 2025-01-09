#!/bin/bash

c=$1

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/select_samples_nuclear.py \
-i /mnt/ceph/users/info/variants/snvs_indels/WES/2022/2022_01_wes_70487_exome/gatk/pvcf/pvcfs_by_chromosome/wes_70487_exome.gatk.chr${c}.vcf.gz \
--overwrite -o chr${c}.ht

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/rekey_nuclear.py \
--i chr${c}.ht \
-o chr${c}_keyed_by_sample.ht



