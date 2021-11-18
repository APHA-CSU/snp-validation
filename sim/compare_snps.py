import pandas as pd
import glob 

from Bio import SeqIO

"""
Calculate performance stats from simulated data
"""

def masked_positions(mask_filepath):
    mask = pd.read_csv(mask_filepath,
        delimiter='\t',
        skiprows=[0,1],
        header=None,
        names=["CHROM", "START", "END", "RPT"]
    )

    # TODO: This is off by one compared to the Emergency Port validation doc
    #    why is that?
    masked_pos = []
    for i, row in mask.iterrows():
        masked_pos.extend(list(range(row['START'], row['END']+1)))

    return masked_pos

def analyse(results_path, sample, mask_filepath):
    """ Compare simulated SNPs data from simuG against btb-seq's snpTable.tab
        If adjust == True: applies the mask to simulated SNPs and pipeline SNPs
        Returns a dictionary of performance stats
    """

    pipeline_directory = glob.glob(results_path + 'btb-seq-results/Results_simulated-reads_*')[0] + '/'
    
    # Load
    simulated_snps = pd.read_csv(results_path + f'simulated-genome/{sample}.simulated.refseq2simseq.map.txt', delimiter='\t')
    pipeline_snps = pd.read_csv(pipeline_directory + f'snpTables/{sample}_snps.tab', delimiter='\t')
    pipeline_genome = load_consensus(pipeline_directory + f'consensus/{sample}_consensus.fas')

    # Trim adapter sequences off pipeline genome
    pipeline_genome = pipeline_genome[29:-31]

    # Extract SNP positions
    simulated_pos = set(simulated_snps['ref_start'].values)
    pipeline_pos = set(pipeline_snps['POS'].values)
    
    # Extract mask positions    
    masked_pos = set(masked_positions(mask_filepath))

    # Extract N's positions in pipeline genome
    n_pos = set([i+30 for i in range(len(pipeline_genome)) if pipeline_genome[i] == 'N'])

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = len(simulated_pos.intersection(pipeline_pos))

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = len(pipeline_pos - simulated_pos)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  len(simulated_pos - pipeline_pos)

    # N - the number of missing sites in the pipeline genome
    n = len(n_pos)

    # The size of the mask
    m = len(masked_pos)

    # TPs in masked regions 
    tp_in_mask = len(masked_pos.intersection(simulated_pos.intersection(pipeline_pos)))

    # FPs in masked regions
    fp_in_mask = len(masked_pos.intersection(pipeline_pos - simulated_pos))

    # FNs in masked regions
    fn_in_mask = len(masked_pos.intersection(simulated_pos - pipeline_pos))

    # Ns in masked regions - should equal the size of the mask
    n_in_mask = len(masked_pos.intersection(n_pos))

    # Compute Performance Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = tp / (tp + fp)

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = tp / (tp + fn)

    # miss rate as FN/(TP + FN)
    miss_rate = fn / (tp + fn)

    # F-score as 2*(precision*recall)/(precision-recall)
    f_score = 2*(precision*sensitivity)/(precision+sensitivity)

    #  total number of errors (FP + FN) per million sequenced bases
    total_errors = fp + fn

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "Ns": n,
        "masked TPs": tp_in_mask,
        "masked FPs": fp_in_mask, 
        "masked FNs": fn_in_mask, 
        "masked Ns": n_in_mask, 
        "mask size": m, # masked Ns and mask size should be equal - could be how the mask positions are calculated.
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "f_score": f_score,
        "total_errors": total_errors,
    }

def load_consensus(path):
    """ Load a consensus file. Returns the first record in a fasta file a string """

    for seq_record in SeqIO.parse(path, "fasta"):
        return str(seq_record.seq)
