import pandas as pd
from Bio import SeqIO

from utils import bcf_summary
import config


"""
    Calculate performance stats from simulated data
"""


def masked_positions(mask_filepath):
    """ Parse mask file path. Returns a list of mask positions """
    
    mask = pd.read_csv(mask_filepath,
        delimiter='\t',
        skiprows=[0,1],
        header=None,
        names=["CHROM", "START", "END", "RPT"]
    )

    masked_pos = []
    for i, row in mask.iterrows():
        # +1 to to adjust for difference in indexing convention between bedfile and snptables
        masked_pos.extend(list(range(row['START']+1, row['END']+1)))

    return masked_pos

def load_consensus(path):
    """ Parse consensus file. Returns the first record in a fasta file a string """

    for seq_record in SeqIO.parse(path, "fasta"):
        return str(seq_record.seq)

def analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path, mask_filepath=None):
    """ Compare simulated SNPs data from simuG against btb-seq's snpTable.tab
        
        Params:
            simulated_snp_path (str): path to simulated snp table
            pipeline_snp_path (str): path to pipeline snp table
            pipeline_genome_path (str): path to pipeline consensus file
            mask_filepath (str): optional path to mask bed file used in pipeline
        
        Returns: 
            (dict) a dictionary of performance stats
    """

    # classify sites as TPs, FPs or FNs
    tp, fp, fn = classify_sites(simulated_snp_path, pipeline_snp_path)

    if mask_filepath:
        # Masked position
        # Load mask positions    
        masked_pos = set(masked_positions(mask_filepath))

        # TPs, FPs and FNs in masked regions 
        fn.update(tp.intersection(masked_pos)) 
        tp = tp.difference(masked_pos)
        fp = fp.difference(masked_pos) 

    # load consensus file
    pipeline_genome = load_consensus(pipeline_genome_path)

    # Extract N's positions in pipeline genome
    # +1 to achieve common indexing with snp tables 
    n_pos = set([i+1 for i in range(len(pipeline_genome)) if pipeline_genome[i] == 'N'])

    # Errors
    n_tp = len(tp)
    n_fp = len(fp)
    n_fn = len(fn)
    n_n_pos = len(n_pos)
    
    # Compute Performance Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = n_tp / (n_tp + n_fp) if (n_tp + n_fp) else float("inf")

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = n_tp / (n_tp + n_fn) if (n_tp + n_fn) else float("inf")

    # miss rate as FN/(TP + FN)
    miss_rate = n_fn / (n_tp + n_fn) if (n_tp + n_fn) else float("inf")

    # F-score as 2*(precision*recall)/(precision-recall)
    f_score = 2*(precision * sensitivity)/(precision + sensitivity)

    # total number of errors (FP + FN) per million sequenced bases
    total_errors = n_fp + n_fn
    
    return {
        "TP": n_tp,
        "FP": n_fp,
        "FN": n_fn,
        "missing_sites": n_n_pos,
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "f_score": f_score,
        "total_errors": total_errors,
    }

def classify_sites(simulated_snp_path, pipeline_snp_path):
    """ Determine TPs, FPs and FNs. TODO: DRY with analyse()"""

    # Load
    simulated_snps = pd.read_csv(simulated_snp_path,  delimiter='\t')
    pipeline_snps = pd.read_csv(pipeline_snp_path, delimiter='\t')

    # Extract SNP positions
    simulated_pos = set(simulated_snps.loc[simulated_snps['variant_type'] == 'SNP', 'ref_start'].values)
    pipeline_pos = set(pipeline_snps.loc[pipeline_snps['TYPE'] == 'SNP', 'POS'].values)
    
    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = simulated_pos.intersection(pipeline_pos)

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = pipeline_pos - simulated_pos

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  simulated_pos - pipeline_pos

    return tp, fp, fn

def site_stats(simulated_snp_path, pipeline_snp_path, bcf_path, mask_filepath):
    """ A data frame that shows stats at each fp/fn site """

    # Classify sites as TP, FP or FN
    tp, fp, fn = classify_sites(simulated_snp_path, pipeline_snp_path)
    
    # Extract mask positions    
    masked_pos = set(masked_positions(mask_filepath))

    # TODO: DRY with analyse
    # FPs and FNs in masked regions 
    fp_in_mask = masked_pos.intersection(fp)
    fn_in_mask = masked_pos.intersection(fn)

    # Summary Bcf
    df = bcf_summary(bcf_path)
    df = df[df.POS.isin(list(fp) + list(fn))]

    # Error Type Column
    df['error_type'] = 'undefined'
    df.loc[df.POS.isin(list(fp)), 'error_type'] = 'fp'
    df.loc[df.POS.isin(list(fn)), 'error_type'] = 'fn'
    
    # AD0 column
    df['AD1/(AD1+AD0)'] = df['AD1'] / (df['AD1'] + df['AD0'])

    # In mask column
    df['in_mask'] = False
    df.loc[df.POS.isin(list(fp_in_mask)), 'in_mask'] = True
    df.loc[df.POS.isin(list(fn_in_mask)), 'in_mask'] = True

    return df

def benchmark(processed_samples, mask_filepath=config.DEFAULT_MASK_PATH):
    """ Assess performance of processed samples. 
        Returns:
            stats_table : sample-wise dataframe that summarises performance
            sitewise_stats: dictionary of DataFrames keyed by sample name
    """
    
    # Initialise
    stats = []
    stats_masked = []
    sitewise_stats = {}

    # Analyse
    for sample in processed_samples:
        simulated_snp_path = sample.genome.snp_table_path
        pipeline_snp_path = sample.sequenced.snp_table_path
        pipeline_genome_path = sample.sequenced.genome_path
        vcf_path = sample.sequenced.vcf_path

        # Performance Stats
        stat = analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path)
        stat["name"] = sample.name
        
        stats.append(stat)

        # Masked Performance Stats
        stat = analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path, mask_filepath)
        stat["name"] = sample.name
        
        stats_masked.append(stat)


        # Site Statistics at fp/fn/tp positions
        site_stat = site_stats(simulated_snp_path, pipeline_snp_path, vcf_path, mask_filepath)

        sitewise_stats[sample.name] = site_stat

    # Combine
    stats_table = pd.DataFrame(stats)
    stats_table_masked = pd.DataFrame(stats_masked)

    return stats_table, stats_table_masked, sitewise_stats
