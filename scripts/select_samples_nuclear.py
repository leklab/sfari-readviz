import argparse
import logging
import hail as hl
from typing import Union
from generate_split_alleles import *

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


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.
    Defaults correspond to gnomAD values.
    """
    return (
        (gq_expr >= adj_gq)
        & hl.cond(gt_expr.is_haploid(), dp_expr >= haploid_adj_dp, dp_expr >= adj_dp)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            .default(
                (ad_expr[gt_expr[0]] / dp_expr >= adj_ab)
                & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            )
        )
    )


def annotate_adj(
    mt: hl.MatrixTable,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid).
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(
        adj=get_adj_expr(
            mt.GT, mt.GQ, mt.DP, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    )



def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Filter genotypes to adj criteria."""
    if "adj" not in list(mt.entry):
        mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    return mt.drop(mt.adj)


def het_hom_hemi_take_expr(mt):
    return hl.struct(S=mt.s, GQ=mt.GQ)


def het_expr(mt):
    return mt.GT.is_het()


def hom_expr(mt):
    return mt.GT.is_diploid() & mt.GT.is_hom_var()

'''
def hemi_expr(mt):
    return hl.or_missing(
        mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar(),
        mt.GT.is_haploid() & (mt.meta.sex == "male") & (mt.GT[0] == 1),
    )
'''

def main(args):

    hl.init(log="./select_samples", master='local[20]',spark_conf={'spark.driver.memory': '8g', 'spark.executor.memory': '8g'})

    '''
    intervals = ['1:100M-200M', '16:29.1M-30.2M', 'X']
    filtered_mt = hl.filter_intervals(
    mt,
    [hl.parse_locus_interval(x, reference_genome='GRCh37') for x in intervals])
    '''
    
    #mt = format_mt(args.input)
    #hl.parse_locus_interval('M:5592-5655',reference_genome='GRCh37')
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('MT:5592-5655',reference_genome='GRCh37')])

    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:5592-5655',reference_genome='GRCh38')])
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:1',reference_genome='GRCh38')])

    #chr14:21,385,610-21,386,182
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrM:500-16000',reference_genome='GRCh38')])

    rg = hl.get_reference('GRCh37')
    grch37_contigs = [x for x in rg.contigs if not x.startswith('GL') and not x.startswith('M')]
    contig_dict = dict(zip(grch37_contigs, ['chr'+x for x in grch37_contigs]))

    mt = hl.import_vcf(args.input,reference_genome='GRCh38',contig_recoding=contig_dict,array_elements_required=False,force_bgz=True,filter='MONOALLELIC')

    #CHD8
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr14:21385610-21386182',reference_genome='GRCh38')])
    #mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr14:21385195-21456127',reference_genome='GRCh38')])

    '''
    exome_intervals = hl.import_locus_intervals('/mnt/home/mlek/ceph/resources/bed_files/exome_evaluation_regions.v1.interval_list', 
                        reference_genome='GRCh38')

    mt = mt.filter_rows(hl.is_defined(exome_intervals[mt.locus]))
    '''

    mt = generate_split_alleles(mt)




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

    mt = mt.select_entries("GT", "GQ", "DP", "AD")

    logger.info("Filtering to entries meeting GQ, DP and other 'adj' thresholds")
    mt = filter_to_adj(mt)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

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

        return -mt.GQ, hl.rand_unif(0, 1, seed=1)


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
