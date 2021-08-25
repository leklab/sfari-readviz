import argparse
import datetime
import hail as hl
import logging
import os

NUM_POSITION_BINS_PER_CHROM = 20000

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


def parse_args():
	"""Parse command line args."""

	p = argparse.ArgumentParser()
	p.add_argument(
		"-i", "--input-ht",
		help="Path of hail table generated by the previous step (select_samples.py)",
		default="gs://gnomad-readviz/v3_and_v3.1/gnomad_v3_readviz_crams.ht",
	)
	p.add_argument(
		"-o", "--output",
		help="Path where the re-keyed hail table will be written. If not specified, output will be written to the same directory as the input",
	)
	p.add_argument(
		"--remove-AC0-variants",
		help="Whether to filter out gnomADv3 AC0 variants. This flag should be used when the input ht contains only gnomADv3 samples",
		type=bool,
	)
	p.add_argument(
		"--variants-ht",
		help="Required if --remove-AC0-variants is used. This is the final released hail table of variants - used to look up which variants ended up being AC0 after QC",
		action="store_true",
		default="gs://gnomad-public-requester-pays/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht",
	)
	p.add_argument(
		"-p", "--output-partitions",
		help="Split data for each sample into this many tsv's. Set to 1 to output a single .tsv.bgz file. "
			 "Setting this >> 1 will produce faster runtimes on large clusters.",
		type=int,
		default=1,
	)
	p.add_argument(
		"-f", "--overwrite-checkpoints",
		action="store_true",
		help="Overwrite/update all hail checkpoints.",
	)
	args = p.parse_args()

	return args


def remove_AC0_variants(ht, variants_ht):
	"""Removes all variants from ht that are AC=0 in the gnomAD public release"""
	gnomad_ht = hl.read_table(variants_ht)
	gnomad_ht_AC0 = gnomad_ht.filter(gnomad_ht.freq.AC[gnomad_ht.globals.freq_index_dict['adj']] == 0, keep=True)

	ht = ht.anti_join(gnomad_ht_AC0)  # remove AC0 variants

	return ht


def explode_table_by_sample(ht):
	"""Explode variant-level ht by sample."""

	logging.info("Input schema:")
	ht.describe()

	ht = ht.annotate(
		samples_w_het_var=ht.samples_w_het_var.map(lambda x: hl.struct(S=x.S, HL=x.HL, het_or_hom_or_hemi=1)),
		samples_w_hom_var=ht.samples_w_hom_var.map(lambda x: hl.struct(S=x.S, HL=x.HL, het_or_hom_or_hemi=2)),
		#samples_w_hemi_var=ht.samples_w_hemi_var.map(lambda x: hl.struct(S=x.S, GQ=x.GQ, het_or_hom_or_hemi=3)),
	)

	#ht = ht.select(samples=ht.samples_w_het_var.extend(ht.samples_w_hom_var.extend(ht.samples_w_hemi_var)))
	ht = ht.select(samples=ht.samples_w_het_var.extend(ht.samples_w_hom_var))

	ht = ht.explode(ht.samples)

	return ht


def rekey_by_sample(ht):
	"""Re-key table by sample id to make subsequent ht.filter(ht.S == sample_id) steps 100x faster"""

	ht = ht.key_by(ht.locus)
	ht = ht.transmute(
		ref=ht.alleles[0],
		alt=ht.alleles[1],
		het_or_hom_or_hemi=ht.samples.het_or_hom_or_hemi,
		#GQ=ht.samples.GQ,
		HL=ht.samples.HL,
		S=ht.samples.S,
	)
	ht = ht.key_by(ht.S)
	ht = ht.transmute(
		chrom=ht.locus.contig.replace("chr", ""),
		pos=ht.locus.position
	)

	logging.info("Schema after re-key by sample:")
	ht.describe()

	return ht


def rekey_by_position_bin(ht):
	ht = ht.key_by(contig=ht.locus.contig, pos_bin=ht.locus.position % NUM_POSITION_BINS_PER_CHROM)
	logging.info("Schema after re-key by position bin:")
	ht.describe()

	return ht


def main():
	args = parse_args()

	#output_file_prefix = os.path.basename(args.input_ht).replace(".ht", "")
	#output_dir = args.output_dir or os.path.dirname(args.input_ht)
	
	#exploded_ht_keyed_by_position_bin = os.path.join(output_dir, f"{output_file_prefix}_keyed_by_position_bin.ht")
	#exploded_ht_checkpoint_path = os.path.join(output_dir, f"{output_file_prefix}_exploded_with_key.ht")
	#exploded_ht_keyed_by_sample_path = os.path.join(output_dir, f"{output_file_prefix}_exploded_keyed_by_sample.ht")
	
	#print(f"Input: {args.input_ht}")
	#print(f"Output dir: {output_dir}")

	ht = hl.read_table(args.input_ht)

	'''
	if args.remove_AC0_variants:
		if not args.variants_ht:
			raise ValueError("--variants-ht must be specified when --remove-AC0-variants flag is used")
		ht = remove_AC0_variants(ht, args.variants_ht)
	'''

	#ht_keyed_by_position_bin = rekey_by_position_bin(ht)

	'''
	ht_keyed_by_position_bin.write(
		exploded_ht_keyed_by_position_bin,
		overwrite=args.overwrite_checkpoints,
	)
	'''

	ht = explode_table_by_sample(ht)

	'''
	ht = ht.checkpoint(
		exploded_ht_checkpoint_path,
		overwrite=args.overwrite_checkpoints,
		_read_if_exists=not args.overwrite_checkpoints,
	)
	'''

	ht_keyed_by_sample = rekey_by_sample(ht)

	ht_keyed_by_sample.write(args.output)

	
	#ht_keyed_by_sample.write(
	#	exploded_ht_keyed_by_sample_path,
	#	overwrite=args.overwrite_checkpoints,
	#)


if __name__ == "__main__":
	main()