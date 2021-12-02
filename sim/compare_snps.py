import pandas as pd
from Bio import SeqIO

from utils import run

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

def load_vcf_tools_output(path):
    pass

def vcf_tools(simulated_vcf_path, pipeline_bcf_path, analysis_file_path):
    # Remove indels from pipeline bcf and decode to vcf
    run(['vcftools', '--bcf', pipeline_bcf_path, '--remove-indels',
           '--recode', '--recode-INFO-all', '--out', analysis_file_path+'pipeline_SNPs'])
    run(['vcftools', '--vcf', simulated_vcf_path, '--diff',
          analysis_file_path+'pipeline_SNPs.recode.vcf', '--diff-site', '--out', analysis_file_path])

def analyse_vcf(simulated_vcf_path, pipeline_vcf_path, analysis_file_path):
    vcf_tools(simulated_vcf_path, pipeline_vcf_path, analysis_file_path)
    load_vcf_tools_output(analysis_file_path)

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