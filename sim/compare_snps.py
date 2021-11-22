import pandas as pd

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

def analyse(simulated_snps, pipeline_snps, mask_filepath):
    """ Compare simulated SNPs data from simuG against btb-seq's snpTable.tab
        If adjust == True: applies the mask to simulated SNPs and pipeline SNPs
        Returns a dictionary of performance stats
    """

    # Load
    simulated = pd.read_csv(simulated_snps, delimiter='\t')
    pipeline = pd.read_csv(pipeline_snps, delimiter='\t')

    # Extract SNP positions
    simulated_pos = set(simulated['ref_start'].values)
    pipeline_pos = set(pipeline['POS'].values)
    masked_pos = set(masked_positions(mask_filepath))

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = len(simulated_pos.intersection(pipeline_pos))

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = len(pipeline_pos - simulated_pos)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  len(simulated_pos - pipeline_pos)

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

    #  total number of errors (FP + FN) per million sequenced bases
    total_errors = fp + fn

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "masked TPs": tp_in_mask,
        "masked FPs": fp_in_mask, 
        "masked FNs": fn_in_mask, 
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "f_score": f_score,
        "total_errors": total_errors
    }

#TODO: This may not be required if we can get away with using bcftools/vcftools
#      for comparisons. Leaving this here for convenience in case those tools aren't suitable 
def load_consensus(path):
    """ Load a consensus file. Returns the first record in a fasta file a string """

    for seq_record in SeqIO.parse(path, "fasta"):
        return str(seq_record)
