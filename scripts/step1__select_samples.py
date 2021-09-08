import argparse
import logging
import hail as hl

#from gnomad.resources import MatrixTableResource
#from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
#from gnomad.utils.filtering import filter_to_adj

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Readviz_prep")
logger.setLevel(logging.INFO)



def het_hom_hemi_take_expr(mt):
    #return hl.struct(S=mt.s, GQ=mt.GQ)
    return hl.struct(S=mt.s, HL=mt.HL)


def het_expr(mt):
    return mt.GT.is_het()


def hom_expr(mt):
    #return mt.GT.is_diploid() & mt.GT.is_hom_var()
    return mt.GT.is_hom_var()

'''
def hemi_expr(mt):
    return hl.or_missing(
        mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar(),
        mt.GT.is_haploid() & (mt.meta.sex == "male") & (mt.GT[0] == 1),
    )
'''


def format_mt(mt_path: str, hl_threshold: float = 0.95, alt_threshold: float = 0.01) -> hl.MatrixTable:
    """Sets HL to zero if below heteroplasmy threshold, adds in genotype based on HL

    :param str mt_path: path to the MatrixTable
    :param float hl_threshold: heteroplasmy level threshold to determine homoplasmic variants
    :param float alt_threshold: heteroplasmy level threshold that must be reached to define the variant as an alternative allele
    :return: MatrixTable with GT
    :rtype: MatrixTable
    """

    # remove alt threshold
    logger.info('Reading in MT...')
    mt = hl.read_matrix_table(mt_path)

    #mt = mt.rename({'VL': 'HL'})

    # TODO: rename artifact-prone-site filter in combine script
    # replace hyphens in filters with underscores
    mt = mt.annotate_rows(filters = mt.filters.map(lambda x: x.replace('-', '_')))

    # convert array of single string to actual array with by splitting string on semicolon
    mt = mt.annotate_entries(FT=hl.str(mt.FT)[2:-2].split(';'))

    # add back in GT based on hl_threshold
    mt = mt.annotate_entries(GT=(hl.case()
        .when((mt.HL < hl_threshold) & (mt.HL > 0.0), hl.parse_call("0/1"))
        .when(mt.HL >= hl_threshold, hl.parse_call("1/1"))
        .when(mt.HL == 0, hl.parse_call("0/0"))
        .default(hl.null(hl.tcall))))

    return mt



def main(args):

    hl.init(log="./select_samples", default_reference="GRCh38")

    '''
    intervals = ['1:100M-200M', '16:29.1M-30.2M', 'X']
    filtered_mt = hl.filter_intervals(
    mt,
    [hl.parse_locus_interval(x, reference_genome='GRCh37') for x in intervals])
    '''
    
    mt = format_mt(args.input)
    #hl.parse_locus_interval('M:5592-5655',reference_genome='GRCh37')
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('MT:5592-5655',reference_genome='GRCh37')])

    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:5592-5655',reference_genome='GRCh38')])
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:1',reference_genome='GRCh38')])

    mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:500-16000',reference_genome='GRCh38')])

    '''
    meta_ht = hl.read_table(args.sample_metadata_ht)
    meta_ht = meta_ht.filter(meta_ht.release & hl.is_defined(meta_ht.project_meta.cram_path))
    meta_ht = meta_ht.select(
        cram_path=meta_ht.project_meta.cram_path,
        crai_path=meta_ht.project_meta.cram_path.replace(".cram", ".cram.crai"),
        sex=meta_ht.project_meta.sex,
    )
    '''

    '''
    mt = MatrixTableResource(args.gnomad_mt).mt()
    mt = hl.MatrixTable(hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True))

    if args.test:
        logger.info("Filtering to chrX PAR1 boundary: chrX:2781477-2781900")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrX:2781477-2781900")])

    meta_join = meta_ht[mt.s]
    mt = mt.annotate_cols(
        meta=hl.struct(
            sex=meta_join.sex,
            cram=meta_join.cram_path,
            crai=meta_join.crai_path,
        )
    )
    logger.info("Filtering to releasable samples with a defined cram path")
    mt = mt.filter_cols(mt.meta.release & hl.is_defined(mt.meta.cram))
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Adjusting samples' sex ploidy")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus,
            mt.GT,
            mt.meta.sex,
            xy_karyotype_str="male",
            xx_karyotype_str="female",
        )
    )
    '''

    #mt = mt.select_entries("GT", "GQ", "DP", "AD")
    mt = mt.select_entries("GT", "HL", "DP")

    '''
    logger.info("Filtering to entries meeting GQ, DP and other 'adj' thresholds")
    mt = filter_to_adj(mt)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    '''

    logger.info(
        f"Taking up to {args.num_samples} samples per site where samples are het, hom_var, or hemi"
    )

    def sample_ordering_expr(mt):
        """It can be problematic for downstream steps when several samples have many times more variants selected
        than in other samples. To avoid this, and distribute variants more evenly across samples,
        add a random number as the secondary sort order. This way, when many samples have an identically high GQ
        (as often happens for common variants), the same few samples don't get selected repeatedly for all common
        variants.
        """

        #return -mt.GQ, hl.rand_unif(0, 1, seed=1)
        return -mt.HL, hl.rand_unif(0, 1, seed=1)


    mt = mt.annotate_rows(
        samples_w_het_var=hl.agg.filter(
            het_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), args.num_samples, ordering=sample_ordering_expr(mt)),
        ),
        samples_w_hom_var=hl.agg.filter(
            hom_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), args.num_samples, ordering=sample_ordering_expr(mt)),
        ),
    )

    ht = mt.rows()

    #ht = ht.select(ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var)
    ht = ht.select(ht.samples_w_het_var, ht.samples_w_hom_var)
    
    ht.write(args.out, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    '''
    parser.add_argument(
        "--test",
        help="Test on chrX", action="store_true",
    )
    '''
    parser.add_argument(
        "--overwrite",
        help="Overwrite if object already exists", action="store_true",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        help="Number of samples to take from each genotype category at each site",
        default=3,
    )
    '''
    parser.add_argument(
        "--gnomad-mt",
        help="Path of the full gnomAD matrix table with genotypes",
        default="gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse_unsplit.repartitioned.mt",
    )
    '''
    '''
    parser.add_argument(
        "--sample-metadata-ht",
        help="Path of the gnomAD sample metadata ht",
        default="gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht",
    )
    '''

    parser.add_argument(
        "--input", "-i",
        help="Path for input matrix table",
        required=True
    )

    parser.add_argument(
        "--out", "-o",
        help="Path for output hail table",
        required=True
    )
    args = parser.parse_args()

    main(args)


"""
Output schema looks like:

----------------------------------------
Global fields:
    None
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'samples_w_het_var': array<struct {
        S: str,
        GQ: int32
    }>
    'samples_w_hom_var': array<struct {
        S: str,
        GQ: int32
    }>
    'samples_w_hemi_var': array<struct {
        S: str,
        GQ: int32
    }>
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------

"""
