import argparse
import datetime
import hail as hl
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def parse_args():
    """Parse command line args."""
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input-ht",
        help="Path of hail table generated by the previous step (rekey.py)",
        default="gs://gnomad-readviz/v3_and_v3.1/gnomad_v3_readviz_crams_exploded_keyed_by_sample.ht",
    )
    p.add_argument(
        "-o", "--output-bucket-path",
        help="Path where the .tsvs for each sample will be written",
        default="gs://gnomad-readviz/v3_and_v3.1/per_sample_tsvs",
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
    p.add_argument(
        "-s", "--start-with-sample-i",
        help="0-based offset into the list of sample ids. This many sample ids from the beginning of the sample ids file will be discarded.",
        type=int,
        default=0,
    )
    p.add_argument(
        "-n", "--n-samples-to-process",
        help="process at most this many sample ids from the sample ids file",
        type=int,
        default=10**9,
    )
    p.add_argument(
        "--sample_ids_path",
        help="A text file containing one sample id per line",
    )
    args = p.parse_args()

    return args


def read_sample_ids(sample_ids_path, start_with_sample_i, n_samples_to_process, n_sample_ids_to_print=10):
    """Read sample ids file.

    Args:
        sample_ids_path (str): sample ids path
        n_sample_ids_to_print (int): log no more than this many sample ids to stdout.
    Return:
        list: sample id strings
    """
    sample_ids = []
    with hl.hadoop_open(sample_ids_path) if sample_ids_path.startswith("gs://" ) else open(sample_ids_path, "rt") as f:
        for i, line in enumerate(f):
            if i < start_with_sample_i:
                continue
            elif i >= start_with_sample_i + n_samples_to_process:
                break

            sample_id = line.rstrip("\n")
            sample_ids.append(sample_id)

            if i <= n_sample_ids_to_print:
                logging.info(sample_id)
                if i == n_sample_ids_to_print and n_sample_ids_to_print > 0:
                    logging.info("...")

    logging.info(f"Parsed {len(sample_ids)} sample ids from {sample_ids_path}")

    return sample_ids


def export_per_sample_tsvs(ht, sample_ids, output_bucket_path, n_partitions_per_sample):
    """Iterate over sample_ids and export all records with that sample id to separate tsv(s).

    Args:
        ht (hail table): output of explode_table_by_sample(..)
        sample_ids (list): collection of sample id strings
        output_bucket_path (str):
        n_partitions_per_sample (int):
    """

    start_time = datetime.datetime.now()
    for i, sample_id in enumerate(sorted(sample_ids)):
        per_sample_ht = ht.filter(ht.S == sample_id, keep=True)

        if i == 0:
            logging.info("Output schema:")
            per_sample_ht.describe()

        if per_sample_ht.count() < 1:
            continue

        if n_partitions_per_sample > 1:
            per_sample_ht = per_sample_ht.naive_coalesce(n_partitions_per_sample)

        # can't remove the key before naive_coalesce for now because of https://github.com/hail-is/hail/issues/8138
        per_sample_ht = per_sample_ht.key_by()

        # re-order columns and drop sample id since it's in the .tsv filename
        per_sample_ht = per_sample_ht.select('chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi', 'GQ')

        if n_partitions_per_sample > 1:
            tsv_output_path = os.path.join(output_bucket_path, sample_id)
            per_sample_ht.export(tsv_output_path, parallel="separate_header")
        else:
            tsv_output_path = os.path.join(output_bucket_path, f"{sample_id}.tsv.bgz")
            per_sample_ht.export(tsv_output_path, header=False)

        logging.info(f"Done exporting {tsv_output_path} to {n_partitions_per_sample} file(s)")

    logging.info(f"Total runtime: {(datetime.datetime.now() - start_time).total_seconds():0.1f} seconds")


def main():
    args = parse_args()

    ht = hl.read_table(args.input_ht)

    sample_ids = read_sample_ids(args.sample_ids_path, args.start_with_sample_i, args.n_samples_to_process)

    export_per_sample_tsvs(ht, sample_ids, args.output_bucket_path, args.output_partitions)


if __name__ == "__main__":
    main()