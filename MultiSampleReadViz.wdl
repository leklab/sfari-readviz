version 1.0

workflow MultiSampleReadVizWorkflow {
	input {

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		#File variants_tsv_bgz = sample[1]
		#File input_bam = sample[2]
		#File input_bai = sample[3]

		Int PADDING_AROUND_VARIANT = 200
		Int SAMPLES_PER_GROUP = 10

		#File inputSamplesFile

		String callset_mt
		File sample_ids
		File sample_bam_tsv

		File perl_script
		File python_script1
		File python_script2
		File python_script3
		File python_script4

	}


	call ReadVizVariantTable {
		input:
			callset_mt = callset_mt,
			python_script = python_script1
	}

	call KeyBySample {
		input:
			readviz_table = ReadVizVariantTable.outfile,
			python_script = python_script2
	}

	call ExportPerSampleTSV {
		input:
			keyBySample_table = KeyBySample.outfile,
			sample_ids = sample_ids,
			python_script = python_script3

	}

	call MakeGATKTSV {
		input:
			samples_with_variants = ExportPerSampleTSV.samples_with_variants,
			sample_bam_tsv = sample_bam_tsv,
			perl_script = perl_script			
	}

	Array[Array[String]] inputSamples = read_tsv(MakeGATKTSV.outfile)

	scatter (sample in inputSamples) {

		#File variants_tsv_bgz = sample[1]
		#File input_bam = sample[2]
		#File input_bai = sample[3]

		String output_prefix = sub(basename(sample[2]), "\\.bam$", "")

		call PrintReadVizIntervals {
			input:
				all_variants_tsv_bgz = ExportPerSampleTSV.outfile,
				variants_tsv_bgz = sample[1],
				input_bam = sample[2],
				input_bai = sample[3],
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_fasta_dict = ref_fasta_dict,
				output_prefix = output_prefix,
				padding_around_variant = PADDING_AROUND_VARIANT,
		}

		call RunHaplotypeCallerBamout {
			input:
				input_bam = PrintReadVizIntervals.output_raw_bam,
				input_bai = PrintReadVizIntervals.output_raw_bai,
				variant_windows_interval_list = PrintReadVizIntervals.variant_windows_interval_list,
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_fasta_dict = ref_fasta_dict,
				output_prefix = output_prefix,
		}

		call ConvertBamToCram {
			input:
				input_bam = RunHaplotypeCallerBamout.output_bamout_bam,
				input_bai = RunHaplotypeCallerBamout.output_bamout_bai,
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_fasta_dict = ref_fasta_dict,
				output_prefix = output_prefix,
		}

		call DeIdentifyBam {
			input:
				sample_id = output_prefix,
				input_cram = ConvertBamToCram.output_bamout_cram,
				input_crai = ConvertBamToCram.output_bamout_cram_crai,
				sample_variants_tsv_bgz = PrintReadVizIntervals.sample_variants_tsv_bgz,
				python_script = python_script4
		}

	}

	Int num_groups = floor(length(DeIdentifyBam.de_identified_bam)/SAMPLES_PER_GROUP) + 1

	scatter (group_num in range(num_groups)) {

		call DetermineBamID {
			input:
				group_num = group_num,
				num_groups = num_groups,
				group_size = SAMPLES_PER_GROUP,
				de_identified_bam = DeIdentifyBam.de_identified_bam
		}		


		call MergeBams {
			input:
				group_num = group_num,
				num_groups = num_groups,
				group_size = SAMPLES_PER_GROUP,
				combined_bamout_id = DetermineBamID.combined_bamout_id,
				de_identified_bam = DeIdentifyBam.de_identified_bam,
				variant_db = DeIdentifyBam.variant_db
		}		
	}


	output {

		File outfile1 = ReadVizVariantTable.outfile
		File outfile2 = KeyBySample.outfile
		File outfile3 = ExportPerSampleTSV.outfile
		File outfile4 = ExportPerSampleTSV.samples_with_variants
		File outfile5 = MakeGATKTSV.outfile

		Array[File] output_hc_bam = RunHaplotypeCallerBamout.output_bamout_bam
		Array[File] output_hc_bai = RunHaplotypeCallerBamout.output_bamout_bai


		Array[File] de_identified_bam = DeIdentifyBam.de_identified_bam
		Array[File] variant_db = DeIdentifyBam.variant_db

		Array[File] gatk_args = MergeBams.gatk_args
		Array[File] combined_bam = MergeBams.combined_bam


		#Array[File] output_raw_bam = PrintReadVizIntervals.output_raw_bam
		#Array[File] output_raw_bai = PrintReadVizIntervals.output_raw_bai

		Array[File] output_bamout_cram = ConvertBamToCram.output_bamout_cram
		Array[File] output_bamout_cram_crai = ConvertBamToCram.output_bamout_cram_crai
		
		#Array[File] output_gvcf = RunHaplotypeCallerBamout.output_gvcf
		#Array[File] output_gvcf_idx = RunHaplotypeCallerBamout.output_gvcf_idx

		# save small intermediate files for debugging
		#Array[File] input_cram_header = PrintReadVizIntervals.input_cram_header
		#Array[File] variant_windows_bed = PrintReadVizIntervals.variant_windows_bed
		#Array[File] variant_windows_interval_list = PrintReadVizIntervals.variant_windows_interval_list

	}
}

task ReadVizVariantTable {

	input {
		String callset_mt
		File python_script
	}

	command {
		conda activate hail_jupyter

		python ${python_script} \
		-i ${callset_mt} \
		--overwrite -o read_viz_mito_variants.ht

		tar cvf read_viz_mito_variants.ht.tar read_viz_mito_variants.ht 

	}

	output {
		File outfile = "read_viz_mito_variants.ht.tar" 
	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}



task KeyBySample {

	input {
		File readviz_table
		File python_script
	}

	command {
		conda activate hail_jupyter

		tar xvf ${readviz_table}

		python ${python_script} \
		--i read_viz_mito_variants.ht \
		-o read_viz_mito_variants_exploded_keyed_by_sample.ht

		tar cvf read_viz_mito_variants_exploded_keyed_by_sample.ht.tar read_viz_mito_variants_exploded_keyed_by_sample.ht

	}

	output {
		File outfile = "read_viz_mito_variants_exploded_keyed_by_sample.ht.tar" 
	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}

task ExportPerSampleTSV {

	input {
		File keyBySample_table
		File sample_ids
		File python_script
	}

	command {
		conda activate hail_jupyter

		tar xvf ${keyBySample_table}

		python ${python_script} \
		-i read_viz_mito_variants_exploded_keyed_by_sample.ht \
		-o ./variants_by_sample \
		--sample_ids_path ${sample_ids}

		tar cvf variants_by_sample.tar variants_by_sample
		ls -1 ./variants_by_sample/ | cut -f1 -d'.' | sort > sample_with_variants.txt
	}

	output {
		File outfile = "variants_by_sample.tar"
		File samples_with_variants = "sample_with_variants.txt" 
	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}


task MakeGATKTSV {

	input {
		File samples_with_variants
		File sample_bam_tsv
		File perl_script			
	}

	command {
		${perl_script} ${sample_bam_tsv} ${samples_with_variants} > sample_subset.tsv
	}

	output {
		File outfile = "sample_subset.tsv"
	}

	runtime {
		cpus: 1
		requested_memory: 4000
	}
}


task DeIdentifyBam {

	input {
		String sample_id
		File input_cram
		File input_crai
		File sample_variants_tsv_bgz
		File python_script
	}

	command {
		conda activate hail_jupyter

		python ${python_script} \
		${sample_id} \
		${input_cram} \
		${sample_variants_tsv_bgz}

		samtools sort -o "${sample_id}.deidentify_output.sorted.bam" "${sample_id}.deidentify_output.bam"
		samtools index "${sample_id}.deidentify_output.sorted.bam"

	}

	output {
		File variant_db = "${sample_id}.deidentify_output.db"
		File de_identified_bam = "${sample_id}.deidentify_output.sorted.bam"
		File de_identified_bai = "${sample_id}.deidentify_output.sorted.bam.bai"		 
	}

	runtime {
		cpus: 4
		requested_memory: 16000
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
	}

	command {
		python << CODE

		import hashlib

		de_identified_bams = ['${sep="','" de_identified_bam}']
		group_bams = de_identified_bams[${group_num}::${num_groups}]

		with open("inputs_${group_num}.list", "w") as fi:
			for i in range(len(group_bams)):
				fi.write(" -I " + group_bams[i]) 

		CODE

		gatk MergeSamFiles \
		--VALIDATION_STRINGENCY SILENT --ASSUME_SORTED --CREATE_INDEX \
		--arguments_file "inputs_${group_num}.list" \
		-O "${combined_bamout_id}.bam"

	}

	output {
		File gatk_args = "inputs_${group_num}.list"
		File combined_bam = "${combined_bamout_id}.bam"

	}

	runtime {
		cpus: 4
		requested_memory: 16000
	}
}


task PrintReadVizIntervals {

	input {
		File all_variants_tsv_bgz
		String variants_tsv_bgz
		File input_bam
		File input_bai

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix

		Int padding_around_variant

	}

	command <<<

		# echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		tar xvf ~{all_variants_tsv_bgz}

		cp ~{variants_tsv_bgz} "~{output_prefix}_variants.tsv.bgz"

		# 1) Convert variants_tsv_bgz to sorted interval list

		zcat ~{variants_tsv_bgz} | awk '{ OFS="\t" } { print("chr"$1, $2 - ~{padding_around_variant} - 1, $2 + ~{padding_around_variant}) }' > variant_windows.bed
		#zcat "${sample} | awk '{ OFS="\t" } { print("chr"$1, $2 - ~{padding_around_variant} - 1, $2 + ~{padding_around_variant}) }' > variant_windows.bed

		# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
		# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
		# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
		# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

		gatk --java-options "-Xms2g" PrintReadsHeader \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			--read-index "~{input_bai}" \
			-O header.bam

		gatk --java-options "-Xms2g" BedToIntervalList \
			--SORT true \
			--SEQUENCE_DICTIONARY header.bam \
			--INPUT variant_windows.bed \
			--OUTPUT variant_windows.interval_list

		# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

		gatk --java-options "-Xms2g" PrintReads \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			--read-index "~{input_bai}" \
			-L variant_windows.interval_list \
			-O "~{output_prefix}.raw.bam"

		ls -lh



		# echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_raw_bam = "${output_prefix}.raw.bam"
		File output_raw_bai = "${output_prefix}.raw.bai"
		File variant_windows_interval_list = "variant_windows.interval_list"

		# save small intermediate files for debugging
		File input_cram_header = "header.bam"
		File variant_windows_bed = "variant_windows.bed"
		File sample_variants_tsv_bgz = "${output_prefix}_variants.tsv.bgz"
	}

	runtime {
		cpus: 1
		requested_memory: 4000
	}
}

task RunHaplotypeCallerBamout {

	input {
		File input_bam
		File input_bai
		File variant_windows_interval_list

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix
	}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		gatk --java-options "-XX:+DisableAttachMechanism -XX:MaxHeapSize=1000m -Xmx7500m" HaplotypeCaller \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			-L ~{variant_windows_interval_list} \
			-bamout "~{output_prefix}.bamout.bam" \
			-O "~{output_prefix}.gvcf"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_bamout_bam = "${output_prefix}.bamout.bam"
		File output_bamout_bai = "${output_prefix}.bamout.bai"
		File output_gvcf = "${output_prefix}.gvcf"
		File output_gvcf_idx = "${output_prefix}.gvcf.idx"
	}

	runtime {
		cpus: 1
		requested_memory: 8000
	}
}

task ConvertBamToCram {

	input {
		File input_bam
		File input_bai

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix

	}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		samtools view -T ~{ref_fasta} -C "~{input_bam}" > "~{output_prefix}.bamout.cram"
		samtools index "~{output_prefix}.bamout.cram" "~{output_prefix}.bamout.cram.crai"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_bamout_cram = "${output_prefix}.bamout.cram"
		File output_bamout_cram_crai = "${output_prefix}.bamout.cram.crai"
	}

	runtime {
		cpus: 1
		requested_memory: 2000
	}
}

