version 1.0

workflow MultiSampleReadVizWorkflow {
	input {

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String chr

		Int PADDING_AROUND_VARIANT
		Int SAMPLES_PER_GROUP

		File python_script4

		File sample_input_file

		String gatk
		String samtools
		String sqlite3

	}

	Array[Array[String]] inputSamples = read_tsv(sample_input_file)

	scatter (sample in inputSamples) {

		String output_prefix = sub(basename(sample[2]), "\\.bam$", "")

		call PerSampleReads {
			input:
				variants_tsv_bgz = sample[1],
				input_bam = sample[2],
				input_bai = sample[3],
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_fasta_dict = ref_fasta_dict,
				output_prefix = output_prefix,
				padding_around_variant = PADDING_AROUND_VARIANT,
				gatk = gatk,
				samtools = samtools,
				python_script = python_script4

		}

	}

	call DetermineNumGroups {
		input:
			num_bam_samples = length(PerSampleReads.de_identified_bam),
			samples_per_group = SAMPLES_PER_GROUP
	}


	scatter (group_num in range(DetermineNumGroups.num_groups)) {

		call DetermineBamID {
			input:
				group_num = group_num,
				num_groups = DetermineNumGroups.num_groups,
				group_size = SAMPLES_PER_GROUP,
				de_identified_bam = PerSampleReads.de_identified_bam
		}		


		call MergeBams {
			input:
				group_num = group_num,
				num_groups = DetermineNumGroups.num_groups,
				group_size = SAMPLES_PER_GROUP,
				combined_bamout_id = DetermineBamID.combined_bamout_id,
				de_identified_bam = PerSampleReads.de_identified_bam,
				variant_db = PerSampleReads.variant_db,
				chr = chr,
				gatk = gatk,
				sqlite3 = sqlite3
		}		
	}

	call combineDBs {
		input:
			chr = chr,
			group_variant_db = MergeBams.combined_db,
			sqlite3 = sqlite3
	}

	output {
		Array[File] combined_bam = MergeBams.combined_bam
		Array[File] combined_bai = MergeBams.combined_bai
		File chr_combined_db = combineDBs.chr_combined_db
	}
}

task DetermineNumGroups {

	input {
		Int num_bam_samples
		Int samples_per_group
	}

	command {
		python << CODE

		import math
		print(math.ceil(${num_bam_samples}/${samples_per_group}))

		CODE
	}

	output {
		Int num_groups = read_int(stdout()) 
	}

	runtime {
		cpus: 1
		requested_memory: 2000
	}
}



task DetermineBamID {

	input {
		Int group_num
		Int num_groups
		Int group_size
		Array[File] de_identified_bam
	}

	command {
		python << CODE

		import hashlib

		de_identified_bams = ['${sep="','" de_identified_bam}']
		group_bams = de_identified_bams[${group_num}::${num_groups}]

		md5_hash = hashlib.md5(", ".join(group_bams).encode('utf-8')).hexdigest()
		combined_bamout_id = "s%d_gs%d_gn%d_gi%04d_h%s" % (len(group_bams),${group_size},${num_groups},${group_num},md5_hash[-9:])

		print(combined_bamout_id)

		CODE

	}

	output {
		String combined_bamout_id = read_string(stdout())
	}

	runtime {
		cpus: 1
		requested_memory: 4000
	}
}


task MergeBams {

	input {
		Int group_num
		Int num_groups
		Int group_size
		String combined_bamout_id
		Array[File] de_identified_bam
		Array[File] variant_db
		String chr
		String gatk
		String sqlite3
	}

	command {
		python << CODE

		def combine_dbs(input_db_paths, set_combined_bamout_id=None, select_chrom=None, sqlite_queries_filename='test.sql', create_index=False):

			sqlite_queries = []
			sqlite_queries.append('CREATE TABLE "variants" ('
				'"id" INTEGER NOT NULL PRIMARY KEY, '
				'"chrom" VARCHAR(2) NOT NULL, '
				'"pos" INTEGER NOT NULL, '
				'"ref" TEXT NOT NULL, '
				'"alt" TEXT NOT NULL, '
				'"zygosity" INTEGER NOT NULL, '
				'"qual" INTEGER NOT NULL, '
				'"combined_bamout_id" TEXT, '
				'"read_group_id" INTEGER NOT NULL);')

			column_names_string = "chrom, pos, ref, alt, zygosity, qual, combined_bamout_id, read_group_id"
			where_clause = 'WHERE chrom="%s"' % select_chrom if select_chrom else ""


			for input_db_path in input_db_paths:
				sqlite_queries.append(
					'ATTACH "%s" as toMerge; ' % input_db_path +
					'BEGIN; ' +
					'INSERT INTO variants (%s) SELECT %s FROM toMerge.variants %s; ' % (column_names_string, column_names_string, where_clause) +
					'COMMIT; ' +
					'DETACH toMerge;')

			if set_combined_bamout_id:
				sqlite_queries.append(
					'UPDATE variants SET combined_bamout_id="%s";' % set_combined_bamout_id)

			if create_index:
				sqlite_queries.append(
					'CREATE INDEX variant_index ON "variants" ("chrom", "pos", "ref", "alt", "zygosity", "qual");')

			sqlite_queries = "\n".join(sqlite_queries)


			with open(sqlite_queries_filename, "wt") as f:
				f.write(sqlite_queries)


		de_identified_bams = ['${sep="','" de_identified_bam}']
		group_bams = de_identified_bams[${group_num}::${num_groups}]

		variant_dbs = ['${sep="','" variant_db}']
		group_dbs = variant_dbs[${group_num}::${num_groups}]


		combine_dbs(input_db_paths=group_dbs, set_combined_bamout_id="chr${chr}_${combined_bamout_id}", sqlite_queries_filename='group${group_num}_queries.sql')

		with open("inputs_${group_num}.list", "w") as fi:
			for i in range(len(group_bams)):
				fi.write(" -I " + group_bams[i]) 

		CODE

		${gatk} MergeSamFiles \
		--VALIDATION_STRINGENCY SILENT --ASSUME_SORTED --CREATE_INDEX \
		--arguments_file "inputs_${group_num}.list" \
		-O "chr${chr}_${combined_bamout_id}.bam"

		${sqlite3} "chr${chr}_${combined_bamout_id}.db" < "group${group_num}_queries.sql"

	}

	output {
		File sql_file = "group${group_num}_queries.sql"
		File gatk_args = "inputs_${group_num}.list"
		File combined_bam = "chr${chr}_${combined_bamout_id}.bam"
		File combined_bai = "chr${chr}_${combined_bamout_id}.bai"
		File combined_db = "chr${chr}_${combined_bamout_id}.db"
	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}

task combineDBs {

	input {
		String chr
		Array[File] group_variant_db
		String sqlite3
	}

	command {
		python << CODE

		def combine_dbs(input_db_paths, set_combined_bamout_id=None, select_chrom=None, sqlite_queries_filename='test.sql', create_index=False):

			sqlite_queries = []
			sqlite_queries.append('CREATE TABLE "variants" ('
				'"id" INTEGER NOT NULL PRIMARY KEY, '
				'"chrom" VARCHAR(2) NOT NULL, '
				'"pos" INTEGER NOT NULL, '
				'"ref" TEXT NOT NULL, '
				'"alt" TEXT NOT NULL, '
				'"zygosity" INTEGER NOT NULL, '
				'"qual" INTEGER NOT NULL, '
				'"combined_bamout_id" TEXT, '
				'"read_group_id" INTEGER NOT NULL);')

			column_names_string = "chrom, pos, ref, alt, zygosity, qual, combined_bamout_id, read_group_id"
			where_clause = 'WHERE chrom="%s"' % select_chrom if select_chrom else ""


			for input_db_path in input_db_paths:
				sqlite_queries.append(
					'ATTACH "%s" as toMerge; ' % input_db_path +
					'BEGIN; ' +
					'INSERT INTO variants (%s) SELECT %s FROM toMerge.variants %s; ' % (column_names_string, column_names_string, where_clause) +
					'COMMIT; ' +
					'DETACH toMerge;')

			if set_combined_bamout_id:
				sqlite_queries.append(
					'UPDATE variants SET combined_bamout_id="%s";' % set_combined_bamout_id)

			if create_index:
				sqlite_queries.append(
					'CREATE INDEX variant_index ON "variants" ("chrom", "pos", "ref", "alt", "zygosity", "qual");')

			sqlite_queries = "\n".join(sqlite_queries)


			with open(sqlite_queries_filename, "wt") as f:
				f.write(sqlite_queries)


		group_variant_dbs = ['${sep="','" group_variant_db}']
		combine_dbs(input_db_paths=group_variant_dbs, sqlite_queries_filename='all_chr${chr}_queries.sql', create_index=True)

		CODE

		${sqlite3} "all_variants.chr${chr}.db" < "all_chr${chr}_queries.sql"

	}

	output {
		File chr_combined_db = "all_variants.chr${chr}.db"
	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}

task PerSampleReads {

	input {
		#File all_variants_tsv_bgz
		String variants_tsv_bgz
		File input_bam
		File input_bai

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix

		Int padding_around_variant
		String gatk

		File python_script
		String samtools

	}

	command <<<

		# echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------


		# 1) Convert variants_tsv_bgz to sorted interval list
		zcat ~{variants_tsv_bgz} | awk '{ OFS="\t" } { print("chr"$1, $2 - ~{padding_around_variant} - 1, $2 + ~{padding_around_variant}) }' > variant_windows.bed

		# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
		# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
		# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
		# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

		~{gatk} --java-options "-Xms2g" PrintReadsHeader \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			--read-index "~{input_bai}" \
			-O header.bam

		~{gatk} --java-options "-Xms2g" BedToIntervalList \
			--SORT true \
			--SEQUENCE_DICTIONARY header.bam \
			--INPUT variant_windows.bed \
			--OUTPUT variant_windows.interval_list

		# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

		~{gatk} --java-options "-Xms2g" PrintReads \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			--read-index "~{input_bai}" \
			-L variant_windows.interval_list \
			-O "~{output_prefix}.raw.bam"

		#task RunHaplotypeCallerBamout
		~{gatk} --java-options "-XX:+DisableAttachMechanism -XX:MaxHeapSize=1000m -Xmx7500m" HaplotypeCaller \
			-R ~{ref_fasta} \
			-I "~{output_prefix}.raw.bam" \
			-L variant_windows.interval_list \
			-bamout "~{output_prefix}.bamout.bam" \
			-O "~{output_prefix}.gvcf"

		#task ConvertBamToCram
		#~{samtools} view -T ~{ref_fasta} -C "~{output_prefix}.bamout.bam" > "~{output_prefix}.bamout.cram"
		#~{samtools} index "~{output_prefix}.bamout.cram" "~{output_prefix}.bamout.cram.crai"

		#task DeIdentifyBam
		python ~{python_script} \
		"~{output_prefix}" \
		"~{output_prefix}.bamout.bam" \
		"~{variants_tsv_bgz}"

		~{samtools} sort -o "~{output_prefix}.deidentify_output.sorted.bam" "~{output_prefix}.deidentify_output.bam"
		#~{samtools} index "~{output_prefix}.deidentify_output.sorted.bam"
		
		# clean up
		rm variant_windows.interval_list
		rm variant_windows.bed

		rm "~{output_prefix}.gvcf"
		rm "~{output_prefix}.gvcf.idx"

		rm header.bam
		rm "~{output_prefix}.raw.bam"
		rm "~{output_prefix}.raw.bai"		
		rm "~{output_prefix}.bamout.bam"
		rm "~{output_prefix}.bamout.bai"
		rm "~{output_prefix}.deidentify_output.bam"

		# echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File variant_db = "~{output_prefix}.deidentify_output.db"
		File de_identified_bam = "~{output_prefix}.deidentify_output.sorted.bam"
		#File de_identified_bai = "~{output_prefix}.deidentify_output.sorted.bam.bai"		 
	}

	runtime {
		cpus: 1
		requested_memory: 4000
	}
}

