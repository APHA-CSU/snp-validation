import pandas as pd
from sklearn.metrics import confusion_matrix

from Bio import SeqIO
# pipeline = pd.read_csv('blah.tab', delimiter='\t')
# simulated = pd.read_csv('AF2122.refseq2simseq.map.txt', delimiter='\t')

def masked_positions(mask_filepath='/home/aaronfishman/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas.rpt.regions'):
    mask = pd.read_csv(mask_filepath, 
        delimiter='\t', 
        skiprows=2,
        header=None,
        names=["name", "start", "end", "mask"]
    )

    # TODO: This is off by one compared to the Emergency Port validation doc
    #    why is that?
    masked_pos = []
    for i, row in mask.iterrows():
        masked_pos.extend(range(row['start'], row['end']+1))
        print('masked_pos', masked_pos[-1])
        break

    return masked_pos

def load_consensus():
    path = '/home/aaronfishman/ebs/pipeline-results/btb-seq-6/pipeline/Results_simulated_reads_22Oct21/consensus/simulated.fas'

    for seq_record in SeqIO.parse(path, "fasta"):
        return seq_record


def analyse(simulated_snps, pipeline_snps):
    simulated = pd.read_csv(simulated_snps, delimiter='\t')
    pipeline = pd.read_csv(pipeline_snps, delimiter='\t')

    simulated_pos = set(simulated['ref_start'].values)
    pipeline_pos = set(pipeline['POS'].values)
    masked_pos = set(masked_positions())
    
    print('before: simulated_pos', len(simulated_pos))
    simulated_pos -= masked_pos
    pipeline_pos -= masked_pos

    print('after: simulated_pos', len(simulated_pos))

    consensus = load_consensus()
    missed = simulated_pos - pipeline_pos

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = len(simulated_pos.intersection(pipeline_pos))

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = len(pipeline_pos - simulated_pos)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  len(simulated_pos - pipeline_pos)

    print("False negatives", simulated_pos - pipeline_pos)


    #### Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = tp / (tp+fp)

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = tp / (tp+fn)

    # miss rate as FN/(TP + FN)
    miss_rate = fn / (tp+fn)

    #  total number of errors (FP + FN) per million sequenced bases
    total_errors = fp+fn

    ### Output
    stats = {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "total_errors": total_errors
    }

    return stats

    a = 1

if __name__ == '__main__':
    root = '/home/aaronfishman/ebs/pipeline-results/btb-seq-6/'

    simulated_snps_filepath = root + 'simulated-genome/simulated.refseq2simseq.map.txt'
    pipeline_snps_filepath = root + 'pipeline/Results_simulated_reads_22Oct21/snpTables/simulated.tab'

    stats = analyse(simulated_snps_filepath, pipeline_snps_filepath)
    print("stats", stats)