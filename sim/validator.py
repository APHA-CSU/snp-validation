import os
import glob
import sys
import argparse
import shutil

import pandas as pd

import compare_snps
from sample import *
from utils import run

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'


def btb_seq(btb_seq_directory, reads_directory, results_directory):
    run(["bash", "./btb-seq", reads_directory,
         results_directory], cwd=btb_seq_directory)

def simulate(
    output_path,
    samples,
    reference_path=DEFAULT_REFERENCE_PATH,
    exist_ok=True
):
    """ Runs a performance test against the pipeline

        Parameters:
            output_path (str): Path to output
            samples (list of Sample objects): Samples for building simulated reads
            reference_path (str): Path to reference genome
            exist_ok (bool): Whether or not to throw an error if a results directory already exists

        Returns:
            None
    """

    # Validate Input
    if os.path.isdir(output_path):
        raise Exception("Output path already exists")
    os.makedirs(output_path)
    
    output_path = os.path.join(output_path, '')

    # Output Directories
    simulated_genome_path = output_path + 'simulated-genome/'
    simulated_reads_path = output_path + 'simulated-reads/'

    # Create Output Directories
    os.makedirs(simulated_genome_path, exist_ok=exist_ok)
    os.makedirs(simulated_reads_path, exist_ok=exist_ok)

    # Simulate Reads
    # TODO:could these function singatures be improved?
    for sample in samples:
        sample.simulate_genome(reference_path, simulated_genome_path + sample.name)
        sample.simulate_reads(simulated_genome_path, simulated_reads_path)


def performance_test(
    output_path,
    btb_seq_path, 
    samples,
    simulated_reads_path=None,  
    exist_ok=True, 
    branch=None,
    light_mode=False
):
    """ Runs a performance test against the pipeline

        Parameters:
            output_path (str): Path to output
            btb_seq_path (str): Path to btb-seq code is stored
            samples (list of Sample objects): Samples on which to run validation on
            simulated_reads_path (str): path to simualted reads
            exist_ok (bool): Whether or not to throw an error if a results directory already exists
            branch (str): Checkout git branch on the repo (default None)

        Returns:
            None
    """

    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    output_path = os.path.join(output_path, '')

    # # Validate Input
    if not os.path.isdir(output_path):
        raise Exception("Output path not found")

    if not os.path.isdir(btb_seq_path):
        raise Exception("Pipeline code repository not found")

    # Output Directories
    results_path = os.path.join(output_path + branch, '')
    os.makedirs(results_path, exist_ok=exist_ok)

    btb_seq_backup_path = results_path + 'btb-seq/'
    btb_seq_results_path = results_path + 'btb-seq-results/'
    site_stats_path = results_path + 'site-stats/'

    # Paths to mask file and simulated reads directory 
    mask_filepath = btb_seq_backup_path + "references/Mycbovis-2122-97_LT708304.fas.rpt.regions"
    if not simulated_reads_path:
        simulated_reads_path = output_path + 'simulated-reads/'

    # TODO: handle dwgsim vcf files. Make sure we are taking into account variants it might generate

    # Create Output Directories
    os.makedirs(btb_seq_results_path, exist_ok=exist_ok)
    os.makedirs(site_stats_path, exist_ok=exist_ok)

    # Backup btb-seq code
    # TODO: exclude the work/ subdirectory from this operation.
    #   This could potentially copy large amounts of data
    #   from the work/ directory nextflow generates
    shutil.copytree(btb_seq_path, btb_seq_backup_path, symlinks=True)

    # TODO: Use nextflow's method of choosing github branches
    #    the method handles automatic pulling of github branches.
    if branch:
        checkout(btb_seq_backup_path, branch)

    # Run the pipeline
    btb_seq(btb_seq_backup_path, simulated_reads_path, btb_seq_results_path)

    # Analyse Results
    # HACK: this could easily break if additional files are present
    pipeline_directory = glob.glob(btb_seq_results_path + 'Results_simulated-reads_*')[0] + '/'

    stats = []
    for sample in samples:
        pipeline_snp_path = pipeline_directory + f'snpTables/{sample.name}_snps.tab'
        if not os.path.exists(pipeline_snp_path):
            pipeline_snp_path = pipeline_directory + f'snpTables/{sample.name}.tab'
        
        if not os.path.exists(pipeline_snp_path):
            raise Exception("Cant Find the pipeline's snps table!!")
        
        simulated_snp_path = output_path + f'simulated-genome/{sample.name}.simulated.refseq2simseq.map.txt'
        pipeline_genome_path = pipeline_directory + f'consensus/{sample.name}_consensus.fas'
        
        if not os.path.exists(pipeline_genome_path):
            pipeline_genome_path = pipeline_directory + f'consensus/{sample.name}.fas'
        
        if not os.path.exists(pipeline_snp_path):
            raise Exception("Cant Find the pipeline's consensus file!!")
        
        # Performance Stats
        stat = compare_snps.analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path, mask_filepath)
        stat["name"] = sample.name
        
        stats.append(stat)

        # Site Statistics at fp/fn/tp positions
        vcf_path = f"{pipeline_directory}/vcf/{sample.name}.vcf.gz"
        site_stats = compare_snps.site_stats(simulated_snp_path, pipeline_snp_path, vcf_path)
        site_stats.to_csv(f"{site_stats_path}/{sample.name}_stats.csv")

    stats_table = pd.DataFrame(stats)

    path = results_path + "stats.csv"
    stats_table.to_csv(path)

    # clean-up
    if light_mode:
        shutil.rmtree(btb_seq_results_path)
        shutil.rmtree(btb_seq_backup_path)

def checkout(repo_path, branch):
    run(["git", "checkout", str(branch)], cwd=repo_path)

def main():
    # Parse
    parser = argparse.ArgumentParser(description="Performance test btb-seq code")
    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("output_path", help="path to performance test results")
    parser.add_argument("--branch", help="name of btb-seq branch to use", default='master')
    parser.add_argument("--ref", "-r", help="optional path to reference fasta", default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--light", "-l", dest='light_mode', help="optional argument to run in light mode", 
                        action='store_true', default=False)
    parser.add_argument("--quick", "-q", help="Run quick samples", action='store_true')

    args = parser.parse_args(sys.argv[1:])

    # Collect Samples
    if args.quick:
        samples = quick_samples()
    else:
        samples = standard_samples()

    # Simulate reads
    simulate(args.output_path, samples)

    # Run
    performance_test(
        args.output_path, 
        args.btb_seq, 
        samples, 
        branch=args.branch,
        light_mode = args.light_mode
    )

if __name__ == '__main__':
    main()