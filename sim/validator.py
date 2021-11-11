import os
import subprocess
import glob
import json
import sys
import argparse
import shutil

import pandas as pd

from compare_snps import analyse
from sample import *
from utils import run

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'


def simulate_reads(
    genome_fasta,
    output_directory,
    sample_name="simulated",
    num_read_pairs=150000,
    read_length=150,
    seed=1,
    outer_distance=330,
    random_dna_probability=0.01,
    rate_of_mutations=0,
    indel_mutation_fraction=0,
    indel_extension_probability=0,
    per_base_error_rate="0" # TODO: default to Ele's reccomendation? 0.001-0.01
):   
    
    # How dwgsim chooses to name it's output fastq files
    output_prefix = output_directory + sample_name
    dwgsim_read_1 = output_prefix + ".bwa.read1.fastq.gz"
    dwgsim_read_2 = output_prefix + ".bwa.read2.fastq.gz"

    # Output fastq filenames compatible with btb-seq
    read_1 = output_prefix + "_S1_R1_X.fastq.gz"
    read_2 = output_prefix + "_S1_R2_X.fastq.gz"

    run([
        "dwgsim",
        "-e", str(per_base_error_rate),
        "-E", str(per_base_error_rate),
        "-i",
        "-d", str(outer_distance),
        "-N", str(num_read_pairs),
        "-1", str(read_length),
        "-2", str(read_length),
        "-r", str(rate_of_mutations),
        "-R", str(indel_mutation_fraction),
        "-X", str(indel_extension_probability),
        "-y", str(random_dna_probability),
        "-H",
        "-z", str(seed),
        genome_fasta,
        output_prefix
    ])

    # Rename output fastq files
    os.rename(dwgsim_read_1, read_1)
    os.rename(dwgsim_read_2, read_2)


def btb_seq(btb_seq_directory, reads_directory, results_directory):
    run(["bash", "./btb-seq", reads_directory,
         results_directory], cwd=btb_seq_directory)


def performance_test(results_path, btb_seq_path, reference_path, exist_ok=True, branch=None):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            results_path (str): Output path to performance test results
            reference_path (str): Path to reference fasta
            exist_ok (bool): Whether or not to throw an error if a results directory already exists
            branch (str): Checkout git branch on the repo (default None)

        Returns:
            None
    """

    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    results_path = os.path.join(results_path, '')

    # # Validate Input
    # if os.path.isdir(results_path):
    #     raise Exception("Output results path already exists")
    # os.makedirs(results_path)

    # if not os.path.isdir(btb_seq_path):
    #     raise Exception("Pipeline code repository not found")

    # Output Directories
    simulated_genome_path = results_path + 'simulated-genome/'
    simulated_reads_path = results_path + 'simulated-reads/'
    btb_seq_backup_path = results_path + 'btb-seq/'
    btb_seq_results_path = results_path + 'btb-seq-results/'

    # Paths to simulated reference genome and simulated SNPs file
    mask_filepath = btb_seq_backup_path + "references/Mycbovis-2122-97_LT708304.fas.rpt.regions"

    # TODO: handle dwgsim vcf files. Make sure we are taking into account variants it might generate

    # Create Output Directories
    os.makedirs(simulated_genome_path, exist_ok=exist_ok)
    os.makedirs(simulated_reads_path, exist_ok=exist_ok)
    os.makedirs(btb_seq_results_path, exist_ok=exist_ok)

    # Backup btb-seq code
    # TODO: exclude the work/ subdirectory from this operation.
    #   This could potentially copy large amounts of data
    #   from the work/ directory nextflow generates
    shutil.copytree(btb_seq_path, btb_seq_backup_path, symlinks=True)

    # TODO: Use nextflow's method of choosing github branches
    #    the method handles automatic pulling of github branches.
    if branch:
        checkout(btb_seq_backup_path, branch)

    # Prepare Genomes
    samples = ["sample1", "sample2"]
    samples = [RandomSample(16000, 1)]
    
    # Simulate Reads
    for sample in samples:
        # TODO: this is ugly, make this human-readable
        sample_name = str(id(sample))

        sample.simulate(reference_path, simulated_genome_path + sample_name + '.')

        # simulate_genome_random_snps(reference_path, simulated_genome_path + sample + '.')

        # # TODO: explicitly path fasta path to simulate
        fasta_path = simulated_genome_path + sample_name + '.simulated.simseq.genome.fa'

        simulate_reads(fasta_path, simulated_reads_path, sample_name=sample_name)

    quit()

    # Run the pipeline
    btb_seq(btb_seq_backup_path, simulated_reads_path, btb_seq_results_path)

    # Analyse Results
    # HACK: this could easily break if additional files are present
    pipeline_directory = glob.glob(btb_seq_results_path + 'Results_simulated-reads_*')[0] + '/'

    stats = []
    for sample in samples:
        simulated_snps = simulated_genome_path + sample + ".simulated.refseq2simseq.map.txt"
        pipeline_snps = pipeline_directory + f'snpTables/{sample}_snps.tab'
        stats.append(analyse(simulated_snps, pipeline_snps, mask_filepath))

    # TODO: add sample name to output

    stats_table = pd.DataFrame(stats)
    stats_table.to_csv(results_path + "stats.csv")

def checkout(repo_path, branch):
    run(["git", "checkout", str(branch)], cwd=repo_path)

def main():
    # Parse
    parser = argparse.ArgumentParser(
        description="Performance test btb-seq code")
    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("results", help="path to performance test results")
    parser.add_argument("--branch", help="name of btb-seq branch to use", default=None)
    parser.add_argument("--ref", "-r", help="optional path to reference fasta", default=DEFAULT_REFERENCE_PATH)

    args = parser.parse_args(sys.argv[1:])

    # Run
    performance_test(args.results, args.btb_seq, args.ref, branch=args.branch)

if __name__ == '__main__':
    main()