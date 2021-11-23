import os
import glob
import sys
import argparse
import shutil

import pandas as pd

from compare_snps import analyse
from sample import *
from utils import run

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'


def btb_seq(btb_seq_directory, reads_directory, results_directory):
    run(["bash", "./btb-seq", reads_directory,
         results_directory], cwd=btb_seq_directory)

# TODO: move these functions into samples

def quick_samples(self):
    return [RandomSample(1)]

def standard_samples(vcf_dir='/mnt/fsx-027/snippy'):
    if not os.path.isdir(vcf_dir):
        raise Exception("Predefined SNP directory not found") 
    
    vcf_dir = os.path.join(vcf_dir, '')

    vcf_filepaths = glob.glob(vcf_dir+'*.vcf')

    samples = []

    samples += [RandomSample(seed=1), RandomSample(seed=666)]

    for filepath in vcf_filepaths:
        samples.append(VcfSample(filepath))    

    return samples

def performance_test(
    results_path, 
    btb_seq_path, 
    reference_path=DEFAULT_REFERENCE_PATH,
    samples=[RandomSample(16000, 1)],
    exist_ok=True, 
    branch=None
):
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
    if os.path.isdir(results_path):
        raise Exception("Output results path already exists")
    os.makedirs(results_path)

    if not os.path.isdir(btb_seq_path):
        raise Exception("Pipeline code repository not found")

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

    # Simulate Reads
    # TODO:could these function singatures be improved?
    for sample in samples:
        sample.simulate_genome(reference_path, simulated_genome_path + sample.name)
        sample.simulate_reads(simulated_genome_path, simulated_reads_path)

    # Run the pipeline
    btb_seq(btb_seq_backup_path, simulated_reads_path, btb_seq_results_path)

    # Analyse Results
    # HACK: this could easily break if additional files are present
    pipeline_directory = glob.glob(btb_seq_results_path + 'Results_simulated-reads_*')[0] + '/'

    stats = []
    for sample in samples:
        pipeline_snps = pipeline_directory + f'snpTables/{sample.name}_snps.tab'

        if not os.path.exists(pipeline_snps):
            pipeline_snps = pipeline_directory + f'snpTables/{sample.name}.tab'

        if not os.path.exists(pipeline_snps):
            raise Exception("Cant Find the pipeline's snps table!!")

        stat = analyse(results_path, sample.name, mask_filepath)
        stat["name"] = sample.name
        
        stats.append(stat)

    stats_table = pd.DataFrame(stats)

    path = results_path + "stats.csv"
    print("***printing to path***")
    stats_table.to_csv(path)

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
    samples = standard_samples()

    performance_test(
        args.results, 
        args.btb_seq, 
        args.ref,
        samples=samples, 
        branch=args.branch
    )

if __name__ == '__main__':
    main()

