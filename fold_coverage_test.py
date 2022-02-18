# Load a 
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import pandas as pd

from samples.simulations.sample import RandomSample, VcfSample
import validator

def run_fold_coverage_test():
    # TODO: remove hardcoded values

    snp_vcf = '/home/aaronfishman/mnt/fsx-027/snippy/AF-21-07727-19.vcf'
    results_path = '/home/aaronfishman/temp-results/fold-coverage-test-5/'
    btb_seq = '/home/aaronfishman/repos/btb-seq'

    genome_len = 4349904
    read_len = 150

    # TODO: DRY up VcfSample contruction

    samples = [
        VcfSample(snp_vcf, num_read_pairs=0.5*10*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*11*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*12*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*13*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*14*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*15*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*16*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*17*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*18*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*19*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*20*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*30*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*40*genome_len/read_len),
        VcfSample(snp_vcf, num_read_pairs=0.5*50*genome_len/read_len),
    ]

    validator.performance_test(results_path, btb_seq, samples=samples)

def plot_result(csv_path='/home/aaronfishman/temp-results/fold-coverage-test-5/stats.csv'):
    # TODO: save plots to output image
    # TODO: remove hardcoded values

    df = pd.read_csv(csv_path)
    fold = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50]

    ax = plt.figure()

    plt.plot(fold, df["FN"])
    plt.plot(fold, df["FN"], 'rx')
    plt.xlabel('Fold Coverage')
    plt.ylabel('False Negatives')
    plt.grid()

    plt.show()
