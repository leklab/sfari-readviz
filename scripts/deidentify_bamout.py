import argparse
import collections
import gzip
import peewee
from pprint import pprint
import logging
import os
import pysam
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

SQLITE_BATCH_SIZE = 10000
FETCH_PADDING_AROUND_VARIANT = 200  # bp


def parse_tsv(tsv_path, progress_bar=False):
    """Parse tsv file, yield: (chrom, pos, ref, alt, zygosity, tsv_path, qual)"""

    with gzip.open(tsv_path, "rt") as f:
        for line in (tqdm(f, unit=f" {tsv_path} lines") if progress_bar else f):
            fields = line.strip().split("\t")
            chrom, pos, ref, alt, zygosity_number, qual = fields
            pos = int(pos)
            zygosity_number = int(zygosity_number)
            #qual = int(qual)
            qual = float(qual)

            yield chrom, pos, ref, alt, zygosity_number, qual


CHROM_ORDER = {k: i for i, k in enumerate(list(map(str, range(1, 23))) + ["X", "Y", "M"])}

def tsv_record_order(tsv_record):
    #tsv_path = item[0]
    chrom, pos, ref, alt, zygosity_number, qual = tsv_record
    return CHROM_ORDER[chrom], pos, len(ref), len(alt), zygosity_number


def interval_union(a, b):
    """Returns a 2-tuple representing the union of two intervals.
    Args:
        a:  2-tuple containing integer coordinates representing a half-open interval.
        b:  2-tuple containing integer coordinates representing the other half-open interval.
    Returns:
        2-tuple containing integer coordinates representing the union(a, b) as a half-open interval.
    """

    return min(a[0], b[0]), max(a[1], b[1])


def do_intervals_intersect(a, b):
    """Returns true if the given 2-tuples overlap.
    Args:
        a:  2-tuple containing integer coordinates representing a half-open interval.
        b:  2-tuple containing integer coordinates representing the other half-open interval.
    Returns:
        True if a and b overlap.
    """

    return a[0] < b[1] and a[1] > b[0]


def read_name_hash(read_name):
    """Takes a read name from an input bam and shortens & obfuscates it"""

    return str(abs(hash(read_name)) % 10**9)  # 9-digit read name


def copy_reads_to_bam(obam, ibam, chrom, pos, ref, alt, zygosity, read_group_id):
    """Copies the reads for the given variant from ibam to obam,
    while discarding extraneous or sensitive information.
    Obfuscates and downsizes read names, discards all tags except an opaque read group
    which is unique for this variant and sample.
    Returns True if the input bam had any reads overlapping this variant.

    Args:
        obam:  output_bam
        ibam:  input bam
        chrom: chromosome (eg. '1' or 'X')
        pos: minrep'ed variant position integer (eg. 12345)
        ref: minrep'ed ref allele (eg. 'A', 'ACT', etc.)
        alt: minrep'ed alt allele (eg. 'GCT', 'C', etc.)
        zygosity: 1 = HET, 2 = HOM, 3 = HEMI
        read_group_id: a hash of the sample id and position
    Return:
        True if the input bam had any reads overlapping this variant
    """

    # compute variant start, end reference coords (half-open)
    variant_start = pos
    variant_end = pos + len(ref)

    # This counter is used as a sanity check that HC added at least one artificial haplotype (typically it adds
    # 2*n of these where n is the number of SNPs in the region).
    artificial_haplotype_counter = 0

    # artificial haplotype coords are half-open (eg. (start=83, end=93) has length 10)
    union_of_artificial_haplotypes_that_overlap_variant = (1e9, 0)  # union of genomic intervals spanned by artificial haplotypes that overlap the variant
    artificial_haplotypes_that_dont_overlap_variant = {}  # maps each artificial haplotype id (eg. HC tag value) to the interval spanned by this artificial haplotype: (r.reference_start, r.reference_end)

    # iterate over the reads
    raw_reads = {}  # maps each artificial haplotype id (eg. HC tag value) to the list of reads assigned to this haplotype (eg. that have this id in their HC tag)
    for r in list(ibam.fetch(f"chr{chrom}", max(1, variant_start - FETCH_PADDING_AROUND_VARIANT), variant_end + FETCH_PADDING_AROUND_VARIANT)):
    #for r in list(ibam.fetch(f"{chrom}", max(1, variant_start - FETCH_PADDING_AROUND_VARIANT), variant_end + FETCH_PADDING_AROUND_VARIANT)):
        tags = dict(r.tags)
        haplotype_id = tags.get('HC', 'reads-without-HC-tag')  # keep reads that aren't asigned to a haplotype
        if tags.get('RG') == "ArtificialHaplotype":
            # handle reads that are actually artificial haplotypes
            artificial_haplotype_counter += 1
            # check whether the artificial haplotype overlaps the variant
            if r.reference_start >= variant_end or r.reference_end <= variant_start:
                # there's no overlap
                artificial_haplotypes_that_dont_overlap_variant[haplotype_id] = (r.reference_start, r.reference_end)
            else:
                union_of_artificial_haplotypes_that_overlap_variant = interval_union(
                    (r.reference_start, r.reference_end),
                    union_of_artificial_haplotypes_that_overlap_variant)

        else:
            # this is a regular read - save it, hashed by the haplotype_id of the haplotype that it was mapped to.
            if haplotype_id not in raw_reads:
                raw_reads[haplotype_id] = []
            raw_reads[haplotype_id].append(r)

    artificial_haplotypes_deleted_counter = 0
    if not raw_reads:
        return False

    # For each artificial haplotype that doesn't overlap the variant, check if it overlaps any of the artificial
    # haplotypes that do overlap the variant. If it does then discard all raw reads that map to it since these reads
    # cause bumps in the coverage plot due to double-counting of the overlapping reads.
    for haplotype_id, artificial_haplotype_that_doesnt_overlap_variant in artificial_haplotypes_that_dont_overlap_variant.items():
        if haplotype_id not in raw_reads:
            continue # skip haplotypes that have no reads mapped to them (this does happen)

        if do_intervals_intersect(
                artificial_haplotype_that_doesnt_overlap_variant,
                union_of_artificial_haplotypes_that_overlap_variant):
            # intersection found, so delete all reads mapping to this haplotype that doesn't overlap the variant
            artificial_haplotypes_deleted_counter += 1
            del raw_reads[haplotype_id]

    # sanity check
    if not raw_reads:
        return False

    read_group_name = f"{chrom}-{pos}-{ref}-{alt}-{zygosity}-{read_group_id}"
    rg_tag = (('RG', read_group_name),)

    found_reads = False
    sorted_reads_iter = sorted((r for hc, reads in raw_reads.items() for r in reads), key=lambda r: r.pos)
    for r in sorted_reads_iter:
        # copy info from r to s
        s = pysam.AlignedSegment()
        s.query_name = read_name_hash(r.query_name)
        s.query_sequence = r.query_sequence
        s.flag = r.flag
        s.reference_id = r.reference_id  # since the bam should only have reads from one chromosome, there will always be just 1 chromosome entry in the header, and so this reference_id can always be 0.
        s.reference_start = r.reference_start
        s.mapping_quality = r.mapping_quality
        s.cigar = r.cigar
        s.next_reference_id = r.next_reference_id
        s.next_reference_start = r.next_reference_start
        s.template_length = r.template_length
        s.query_qualities = r.query_qualities
        s.tags = rg_tag

        obam.write(s)
        found_reads = True

    return found_reads


def generate_deidentified_bam(sample_id, input_bam_path, tsv_record_iterator, tsv_records_to_skip=None):
    """Reads variants from the given tsv_record_iterator and creates the deidentified bam
     by copying reads from the input bam while removing any extraneous or sensitive info
     from the reads and bam header.

    Args:
        sample_id: sample id
        input_bam_path: input bam path
        tsv_record_iterator: returns parsed tsv records
        tsv_records_to_skip: optional Set of records to exclude from the deidentified bam
    """
    try:
        ibam = pysam.AlignmentFile(input_bam_path, "rb")
    except Exception as e:
        raise Exception(f"Unable to open: {input_bam_path}: {e}")

    # create output bam header
    reference_sequences = []
    ibam_header_dict = dict(ibam.header)  # read entire header dict into memory here to avoid multiple i/o operations
    for reference_id in range(len(ibam_header_dict['SQ'])):
        d = {}
        reference_sequences.append(d)
        for key, value in ibam_header_dict['SQ'][reference_id].items():
            if key in ["SN", "LN"]:
                d[key] = value

    header = {
        'HD': {'VN': '1.4', 'SO': 'coordinate'},
        'SQ': reference_sequences,
        'RG': [],
    }

    # open output bam file
    bam_output_path = f"{sample_id}.deidentify_output.bam"
    with pysam.AlignmentFile(bam_output_path, "wb", header=header) as obam:

        counter = 0
        prev_chrom = None
        for record in tsv_record_iterator:
            chrom, pos, ref, alt, zygosity_number, qual = record

            if tsv_records_to_skip and record in tsv_records_to_skip:
                logging.info(" ".join(map(str, ["Skipping", chrom, pos, ref, alt, zygosity_number, qual])) +
                    f" since it's in the list of {len(tsv_records_to_skip)} records to skip")
                continue

            if chrom != prev_chrom:
                logging.info(f"Processing chrom: {chrom} from {record}")
                prev_chrom = chrom

            # security: python hash is essentially impossible to invert, so there's no danger of someone
            # recovering the sample id by inverting the hash function. Since python3.3, the hash function
            # is salted with a random seed, which prevents dictionary/rainbow attacks.
            # Apply % 10**6 to save space.
            read_group_id = abs(hash((chrom, pos, ref, alt, zygosity_number, sample_id))) % 10**6

            # fetch reads and write to combined file
            succeeded = copy_reads_to_bam(
                obam, ibam, chrom, pos, ref, alt, zygosity_number, read_group_id)

            if succeeded:
                counter += 1
                yield chrom, pos, ref, alt, zygosity_number, qual, read_group_id
            else:
                logging.warning(" ".join(map(str, ["Skipping", chrom, pos, ref, alt, zygosity_number, qual,
                    "because failed to copy any reads from bam"])))

    logging.info(f"Added reads for {counter} records")

    #os.system(f"samtools sort -o {bamout_group_name}.sorted.bam {bam_output_path}")
    #os.system(f"samtools index {bamout_group_name}.sorted.bam")


def write_to_database(sample_id, db_record_iterator):
    """Generates a sqlite database that allows queries for the number of available tracks for a
    given variant, as well as what bam files contain those reads.
    Supports queries like:

       select chrom, pos, ref, alt, count(*) as c, group_concat(combined_bamout_id), group_concat(qual) from (
           select * from variants order by (qual, read_group_id) desc
       ) group by chrom, pos, ref, alt
    """

    # create sqlite database
    db_path = f"{sample_id}.deidentify_output.db"
    if os.path.isfile(db_path):
        os.remove(db_path)
    sqlite_db = peewee.SqliteDatabase(db_path, autocommit=False)
    class Variants(peewee.Model):
        chrom = peewee.CharField(max_length=2, null=False)
        pos = peewee.IntegerField(null=False)
        ref = peewee.TextField(null=False)
        alt = peewee.TextField(null=False)
        zygosity = peewee.IntegerField(null=False)
        qual = peewee.IntegerField(null=False)
        combined_bamout_id = peewee.TextField(null=True)
        read_group_id = peewee.IntegerField(null=False)

        class Meta:
            database = sqlite_db
            indexes = (
                (('chrom', 'pos', 'ref', 'alt', 'zygosity', 'qual'), False),  # True means unique index
            )

    Variants.create_table(fail_silently=False)

    counter = 0
    db_batch = []
    for record in db_record_iterator:
        chrom, pos, ref, alt, zygosity_number, qual, read_group_id = record
        db_batch.append(Variants(**{
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "zygosity": zygosity_number,
            "qual": qual,
            "combined_bamout_id": None,
            "read_group_id": read_group_id,
        }))

        if len(db_batch) > SQLITE_BATCH_SIZE:
            counter += len(db_batch)
            Variants.bulk_create(db_batch)
            db_batch = []

    Variants.bulk_create(db_batch)
    counter += len(db_batch)

    logging.info(f"Wrote {counter} records to {db_path}")


def main(args):

    tsv_records_to_skip = None

    #if args.input_tsv_path_of_variants_to_skip:
        #tsv_records_to_skip = set(parse_tsv(args.input_tsv_path_of_variants_to_skip, progress_bar=True))

    # assemble pipeline
    tsv_record_iterator = sorted(parse_tsv(args.input_tsv_path, progress_bar=False), key=tsv_record_order)

    combined_bam_generator = generate_deidentified_bam(
        args.sample_id, args.input_bam_path, tsv_record_iterator, tsv_records_to_skip=tsv_records_to_skip)
    
    write_to_database(args.sample_id, combined_bam_generator)


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("-x", "--input-tsv-path-of-variants-to-skip", help="(optional) if specified, variants in this tsv will be excluded from the output bam")
    p.add_argument("sample_id", help="sample id")
    p.add_argument("input_bam_path", help="path of bamout to deidentify")
    p.add_argument("input_tsv_path", help="tsv file for the input bam as generated by step3__export_per_sample_tsvs.py")

    args = p.parse_args()

    #if args.input_tsv_path_of_variants_to_skip and not os.path.isfile(args.input_tsv_path_of_variants_to_skip):
    #    p.error(f"path {args.args.input_tsv_path_of_variants_to_skip} in {args.args.input_tsv_path_of_variants_to_skip} doesn't exist")

    if not os.path.isfile(args.input_tsv_path):
        p.error(f"path {args.input_tsv_path} in {args.input_tsv_path} doesn't exist")

    if not os.path.isfile(args.input_bam_path):
        p.error(f"path {args.input_bam_path} in {args.input_bam_path} doesn't exist")


    main(args)
