from re import I
from numpy.lib import utils
import pandas as pd
import glob 
import os
from Bio import SeqIO

from utils import bcf_summary

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

def analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path, mask_filepath):
    """ Compare simulated SNPs data from simuG against btb-seq's snpTable.tab
        
        Params:
            simulated_snp_path (str): path to simulated snp table
            pipeline_snp_path (str): path to pipeline snp table
            pipeline_genome_path (str): path to pipeline consensus file
            mask_filepath (str): path to mask bed file used in pipeline
        
        Returns: 
            (dict) a dictionary of performance stats
    """
    # Load
    simulated_snps = pd.read_csv(simulated_snp_path,  delimiter='\t')
    pipeline_snps = pd.read_csv(pipeline_snp_path, delimiter='\t')
    pipeline_genome = load_consensus(pipeline_genome_path)

    # Extract SNP positions
    simulated_pos = set(simulated_snps.loc[simulated_snps['variant_type'] == 'SNP', 'ref_start'].values)
    pipeline_pos = set(pipeline_snps.loc[pipeline_snps['TYPE'] == 'SNP', 'POS'].values)
    
    # Extract mask positions    
    masked_pos = set(masked_positions(mask_filepath))

    # Extract N's positions in pipeline genome
    # +1 to achieve common indexing with snp tables 
    n_pos = set([i+1 for i in range(len(pipeline_genome)) if pipeline_genome[i] == 'N'])

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = len(simulated_pos.intersection(pipeline_pos))

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = len(pipeline_pos - simulated_pos)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  len(simulated_pos - pipeline_pos)

    # N - the number of missing sites in the pipeline genome
    n = len(n_pos)

    # TPs in masked regions 
    tp_in_mask = len(masked_pos.intersection(simulated_pos.intersection(pipeline_pos)))

    # FPs in masked regions
    fp_in_mask = len(masked_pos.intersection(pipeline_pos - simulated_pos))

    # FNs in masked regions
    fn_in_mask = len(masked_pos.intersection(simulated_pos - pipeline_pos))

    # Compute Performance Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = tp / (tp + fp) if (tp + fp) else float("inf")

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = tp / (tp + fn) if (tp + fn) else float("inf")

    # miss rate as FN/(TP + FN)
    miss_rate = fn / (tp + fn) if (tp + fn) else float("inf")

    # F-score as 2*(precision*recall)/(precision-recall)
    f_score = 2*(precision*sensitivity)/(precision+sensitivity)

    # total number of errors (FP + FN) per million sequenced bases
    total_errors = fp + fn

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "misssing_sites": n,
        "TPs_in_mask": tp_in_mask,
        "FPs_in_mask": fp_in_mask, 
        "FNs_in_mask": fn_in_mask, 
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "f_score": f_score,
        "total_errors": total_errors
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

def site_stats(simulated_snp_path, pipeline_snp_path, bcf_path):
    print("simulate", "*******")
    print("simulate", simulated_snp_path)
    print("pipeline_snp_path", pipeline_snp_path)
    print("bcf_path", bcf_path)

    tp, fp, fn = classify_sites(simulated_snp_path, pipeline_snp_path)

    # Summary Bcf
    df = bcf_summary(bcf_path)
    df = df[(df.POS.isin(list(fp) + list(fn)))]   

    # Error Type Column
    df['error_type'] = 'undefined'
    df.loc[df.POS.isin(list(fp)), 'error_type'] = 'fp'
    df.loc[df.POS.isin(list(fn)), 'error_type'] = 'fn'
    
    # AD0 column
    df['AD1/(AD1+AD0)'] = df['AD1'] / (df['AD1'] + df['AD0'])

    return df

if __name__ == '__main__':
    pipeline_snps = '/home/aaronfishman/temp/sites3/btb-seq-results/Results_simulated-reads_08Dec21/snpTables/RandomSample-snps16000-indels1600-seed1_snps.tab'
    simulated_snps = '/home/aaronfishman/temp/sites3/simulated-genome/RandomSample-snps16000-indels1600-seed1.simulated.refseq2simseq.map.txt'
    # bcf_path = '/home/aaronfishman/temp/ad0-not/btb-seq-results/Results_simulated-reads_07Dec21/filteredBcf/VcfSample-AF-12-00945-19-seed1-error0.001-0.01-readpairs-144997_filtered.bcf'
    bcf_path = '/home/aaronfishman/temp/sites3/btb-seq-results/Results_simulated-reads_08Dec21//vcf/RandomSample-snps16000-indels1600-seed1.vcf.gz'
    print(site_stats(simulated_snps, pipeline_snps, bcf_path))
