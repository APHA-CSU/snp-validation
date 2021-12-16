import os
import glob
import sys
import argparse
import shutil

import pandas as pd

import sample_sets
import utils
import compare_snps
import sequenced
import processed

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'
DEFAULT_MASK_PATH = './Mycbovis-2122-97_LT708304.fas.rpt.regions'

def simulate(
    samples,
    genomes_path,
    reads_path, 
):
    """ Simulate a number of samples

        Parameters:
            output_path (str): Path to output
            samples (list of Sample objects): Samples for building simulated reads
            reference_path (str): Path to reference genome
            exist_ok (bool): Whether or not to throw an error if a results directory already exists

        Returns:
            simulate_reads_path (str): Path to simulated reads 
    """

    # Validate Input
    reads_path = os.path.join(reads_path, '')
    genomes_path = os.path.join(genomes_path, '')
    os.makedirs(genomes_path, exist_ok=True)
    os.makedirs(reads_path, exist_ok=True)

    # Simulate genome and reads
    simulated_samples = []
    for sample in samples:
        simulated = sample.simulate_genome(genomes_path + sample.name)
        sample.simulate_reads(simulated.genome_path, reads_path)
        simulated_samples.append(simulated)

    # TODO: pass reads path into function rather than generating inside
    return simulated_samples

def sequence(btb_seq_path, reads_path, results_path):
    """ Sequence reads with btb-seq """
    # Validate
    if not os.path.exists(btb_seq_path):
        raise Exception("btb-seq code does not exist: ", btb_seq_path)

    if not os.path.exists(reads_path):
        raise Exception("reads path does not exist: ", reads_path)

    os.makedirs(results_path, exist_ok=True)

    # Sequence
    utils.run(["bash", "./btb-seq", reads_path, results_path], cwd=btb_seq_path)

    # Result directory
    # TODO: handle when glob does not return a unique path
    return glob.glob(results_path + '/Results_*')[0] + '/'

def performance_test(
    btb_seq_path,
    output_path,
    samples,
    light_mode=False
):
    """ Runs a performance test against the pipeline

        Parameters:
            output_path (str): Path to output - for simulations
            btb_seq_path (str): Path to btb-seq code is stored
            samples (list of Sample objects): Samples on which to run validation on
            results_path (str): path to results - if None defaults to output_path
            simulated_reads_path (str): path to simualted reads
            exist_ok (bool): Whether or not to throw an error if a results directory already exists
            branch (str): Checkout git branch on the repo (default None)
            light_mode (bool): switch to delete non-essential analysis files (default False)

        Returns:
            None
    """
    # Define output paths
    genomes_path = os.path.join(output_path, 'genomes')
    reads_path = os.path.join(output_path, 'reads')
    btb_seq_backup_path = os.path.join(output_path, 'btb-seq')
    results_path = os.path.join(output_path, 'sequenced')
    stats_path = os.path.join(output_path, 'stats')

    # Initialise
    # TODO: exclude the work/ subdirectory from this operation to save space
    shutil.copytree(btb_seq_path, btb_seq_backup_path, symlinks=True)
    os.makedirs(stats_path, exist_ok=True)

    # Simulate
    simulated_samples = simulate(samples, genomes_path, reads_path)
    results_path = sequence(btb_seq_path, reads_path, results_path)
    sequenced_samples = sequenced.from_results_dir(results_path)
    processed_samples = processed.from_list(simulated_samples, sequenced_samples)

    stats, site_stats = compare_snps.benchmark(processed_samples)

    # Save
    stats.to_csv(output_path + '/stats.csv')

    for name, df in site_stats.items():
        df.to_csv(stats_path + f'/{name}_stats.csv')

def main():
    # Parse
    parser = argparse.ArgumentParser(description="Performance test btb-seq code")
    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("output_path", help="path to performance test results")
    parser.add_argument("--branch", help="name of btb-seq branch to use", default=None)
    parser.add_argument("--ref", "-r", help="optional path to reference fasta", default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--light", "-l", 
        dest='light_mode', 
        help="optional argument to run in light mode", 
        action='store_true', default=False
    )
    parser.add_argument("--quick", "-q", help="Run quick samples", action='store_true')

    args = parser.parse_args(sys.argv[1:])

    # Collect Samples
    if args.quick:
        samples = sample_sets.quick_samples()

    else:
        samples = sample_sets.standard_samples()

    # Set branch
    if args.branch:
        utils.checkout(args.btb_seq, args.branch)

    # Simulate reads
    simulated_reads_path = simulate(args.output_path, samples)

    # Run
    performance_test(
        args.output_path, 
        args.btb_seq, 
        samples, 
        simulated_read_path=simulated_reads_path,
        light_mode = args.light_mode
    )

if __name__ == '__main__':
    btb_seq_path = '/home/aaronfishman/repos/btb-seq/'
    output_path = '/home/aaronfishman/temp/cleanup/'

    samples = [sample_sets.standard_samples()[0], sample_sets.RandomSample('/home/aaronfishman/tinygenome.fas', num_snps=0, num_indels=0)]

    performance_test(
        btb_seq_path,
        output_path,
        samples,
        light_mode=False
    )
    quit()


    ###### main()
    samples = [standard_samples()[0], RandomSample('/home/aaronfishman/tinygenome.fas', num_snps=0, num_indels=0)]
        
    genomes_path = '/home/aaronfishman/temp/genomes/'
    reads_path = '/home/aaronfishman/temp/reads/'
    results_path = '/home/aaronfishman/temp/results/'
    btb_seq_path = '/home/aaronfishman/repos/btb-seq/'

    simulated_samples = simulate(samples, genomes_path, reads_path)
    results_path = sequence(btb_seq_path, reads_path, results_path)
    sequenced_samples = sequenced.from_results_dir(results_path)
    processed_samples = processed.from_list(simulated_samples, sequenced_samples)

    stats, site_stats = benchmark(processed_samples)


    a = 1
