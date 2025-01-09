#!/bin/bash

c=$1

#below should be a loop but don't have time to figure out and test in bash

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_1_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 0 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_2_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 10000 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_3_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 20000 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_4_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 30000 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_5_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 40000 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_6_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 50000 \
-n 10000

python /mnt/home/mlek/ceph/dev/sfari-readviz/scripts/export_per_sample_tsvs_nuclear.py \
-i chr${c}_keyed_by_sample.ht \
-o /mnt/home/mlek/mlek/readviz/iWES/chr${c}_7_tsvs \
--sample_ids_path iWES_sample_ids.txt \
-s 60000 \
-n 10488
