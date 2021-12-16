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

import config

def simulate(
    samples,
    genomes_path,
    reads_path, 
):
    """ Simulate a number of samples

        Parameters:
            samples (list of Sample objects): Samples for generating simulated reads
            genomes_path (string): path to directory containing output genomes
            reads_path (string): path to directory containing output fastq reads

        Returns:
            simulated_samples (str): Path to simulated reads 
    """

    # Initialise
    reads_path = os.path.join(reads_path, '')
    genomes_path = os.path.join(genomes_path, '')
    os.makedirs(genomes_path, exist_ok=False)
    os.makedirs(reads_path, exist_ok=False)

    # Simulate genome and reads
    simulated_samples = []
    for sample in samples:
        simulated = sample.simulate_genome(genomes_path + sample.name)
        sample.simulate_reads(simulated.genome_path, reads_path)
        simulated_samples.append(simulated)

    return simulated_samples

def sequence(btb_seq_path, reads_path, results_path):
    """ Sequence reads with btb-seq """
    # Validate
    if not os.path.exists(btb_seq_path):
        raise Exception("btb-seq code does not exist: ", btb_seq_path)

    if not os.path.exists(reads_path):
        raise Exception("reads path does not exist: ", reads_path)

    os.makedirs(results_path, exist_ok=False)

    # Sequence
    utils.run(["bash", "./btb-seq", reads_path, results_path], cwd=btb_seq_path)

    # Result directory
    # TODO: handle when glob does not return a unique path
    return glob.glob(results_path + '/Results_*')[0] + '/'

def pipeline(
    btb_seq_path,
    output_path,
    samples,
    light_mode=False
):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            output_path (str): Path to output - for simulations
            samples (list of Sample objects): Samples on which to run validation on
            light_mode (bool): keep only output statistics to save disk space (default False)

        Returns:
            None
    """
    # Define output paths
    genomes_path = os.path.join(output_path, 'genomes')
    reads_path = os.path.join(output_path, 'reads')
    btb_seq_backup_path = os.path.join(output_path, 'btb-seq')
    results_path_parent = os.path.join(output_path, 'sequenced')
    stats_path = os.path.join(output_path, 'stats')

    # Initialise
    # TODO: exclude the work/ subdirectory from this operation to save space
    os.makedirs(output_path, exist_ok=False)
    shutil.copytree(btb_seq_path, btb_seq_backup_path, symlinks=True)
    os.makedirs(stats_path, exist_ok=False)

    # Simulate
    simulated_samples = simulate(samples, genomes_path, reads_path)
    results_path = sequence(btb_seq_path, reads_path, results_path_parent)
    sequenced_samples = sequenced.from_results_dir(results_path)
    processed_samples = processed.from_list(simulated_samples, sequenced_samples)

    stats, site_stats = compare_snps.benchmark(processed_samples)

    # Save
    stats.to_csv(output_path + '/stats.csv')

    for name, df in site_stats.items():
        df.to_csv(stats_path + f'/{name}_stats.csv')

    # Cleanup
    if light_mode:
        shutil.rmtree(genomes_path)
        shutil.rmtree(reads_path)
        shutil.rmtree(btb_seq_backup_path)
        shutil.rmtree(results_path)

def main():
    parser = argparse.ArgumentParser(prog="Performance test btb-seq code")
    subparsers = parser.add_subparsers(help='sub-command help')

    # Full SNP Validation pipeline
    subparser = subparsers.add_parser('pipeline', help='Run the full validation pipeline: simulation, sequencing and benchmarking')
    subparser.add_argument("btb_seq_path", help="path to btb-seq code")
    subparser.add_argument("output_path", help="path to performance test results")
    subparser.add_argument("--light", "-l", action='store_true', dest='light_mode', help="optional argument to run in light mode")
    subparser.add_argument("--quick", "-q", help="Run quick samples", action='store_true')
    subparser.set_defaults(func=pipeline)

    # Simulation
    subparser = subparsers.add_parser('simulate', help='Simulate genomes and reads')
    subparser.add_argument("genomes_path", help="path to directory containing output genomes")
    subparser.add_argument("reads_path", help="path to directory containing output reads")
    subparser.add_argument("--quick", "-q", help="Run quick samples", action='store_true')
    subparser.set_defaults(func=simulate)

    # TODO: benchmarking. Might be useful to be able to reprocess analysis if it changes

    # Parse
    kwargs = vars(parser.parse_args())

    if not kwargs:
        parser.print_help()
        return

    if "quick" in kwargs:
        kwargs["samples"] = sample_sets.quick_samples() if kwargs["quick"] else sample_sets.standard_samples()
        del kwargs["quick"]

    # Run chosen option
    func = kwargs.pop("func")
    func(**kwargs)

if __name__ == '__main__':
    main()
