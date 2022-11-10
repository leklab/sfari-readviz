version 1.0

import "MultiSampleReadViz.wdl" as single_wf

workflow MultiChromReadVizWorkflow {
	input {

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		Array[String] chr_list

		Int PADDING_AROUND_VARIANT = 200
		Int SAMPLES_PER_GROUP = 500

		File python_script4

		Array[File] sample_input_files

		String gatk
		String samtools
		String sqlite3

	}

	scatter (i in range(length(sample_input_files))) {
		call single_wf.MultiSampleReadVizWorkflow {
			input:
				ref_fasta = ref_fasta,
				ref_fasta_fai = ref_fasta_fai,
				ref_fasta_dict = ref_fasta_dict,
				PADDING_AROUND_VARIANT = PADDING_AROUND_VARIANT,
				SAMPLES_PER_GROUP = SAMPLES_PER_GROUP,
				python_script4 = python_script4,
				gatk = gatk,
				samtools = samtools,
				sqlite3 = sqlite3,
				sample_input_file = sample_input_files[i],
				chr = chr_list[i]
		}
	}

	output {

		Array[Array[File]] combined_bam = MultiSampleReadVizWorkflow.combined_bam
		Array[Array[File]] combined_bai = MultiSampleReadVizWorkflow.combined_bai
		Array[File] chr_combined_db = MultiSampleReadVizWorkflow.chr_combined_db
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







