import glob
import os

import pandas as pd

from sample import NamedSample
from compare_snps import analyse

root = '/home/aaronfishman/temp-results/ofat-25/'
root = '/home/aaronfishman/mnt/fsx-027/ofat-25/'

branch_paths = sorted(glob.glob(root + './*'))

for branch_path in branch_paths:
    branch_name = os.path.basename(branch_path)

    # Sample names
    raw_sample_paths = glob.glob(branch_path + '/simulated-reads/*_S1_R2_X.fastq.gz')
    raw_sample_names = [os.path.basename(p).split('.')[0][:-8] for p in raw_sample_paths]

    sample_names = list(set(raw_sample_names))


    # Performa analysis

    # Output Directories
    results_path = branch_path + '/'
    simulated_genome_path = results_path + 'simulated-genome/'
    simulated_reads_path = results_path + 'simulated-reads/'
    btb_seq_backup_path = results_path + 'btb-seq/'
    btb_seq_results_path = results_path + 'btb-seq-results/'
    mask_filepath = btb_seq_backup_path + "references/Mycbovis-2122-97_LT708304.fas.rpt.regions"

    samples = []
    for name in sample_names:
        sample = NamedSample(name)
        samples.append(sample)

    samples = sorted(samples, key=lambda x: x.name)

    # Analyse Results
    # HACK: this could easily break if additional files are present
    pipeline_directory = glob.glob(btb_seq_results_path + '/Results_simulated-reads_*')[0] + '/'

    stats = []
    for sample in samples:
        simulated_snps = simulated_genome_path + sample.name + ".simulated.refseq2simseq.map.txt"
        
        pipeline_snps = pipeline_directory + f'snpTables/{sample.name}_snps.tab'

        if not os.path.exists(pipeline_snps):
            pipeline_snps = pipeline_directory + f'snpTables/{sample.name}.tab'

        if not os.path.exists(pipeline_snps):
            raise Exception("Cant Find the pipeline's snps table!!")

        stat = analyse(simulated_snps, pipeline_snps, mask_filepath)
        stat["name"] = sample.name
        
        stats.append(stat)

    stats_table = pd.DataFrame(stats)
    cols = stats_table.columns.tolist()
    stats_table = stats_table[cols[-1:] + cols[:-1]]

    stats_table.to_csv(results_path + '/stats.csv')

# print(branches)