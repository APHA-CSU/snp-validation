import os
import glob
import sys
import argparse
import shutil

import pandas as pd

import compare_snps
from sample import *
import utils
import processed_sample
import sequenced_sample
import simulated_genome

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'
DEFAULT_MASK_PATH = './Mycbovis-2122-97_LT708304.fas.rpt.regions'

def btb_seq(btb_seq_directory, reads_directory, results_directory):
    utils.run(["bash", "./btb-seq", reads_directory, results_directory], cwd=btb_seq_directory)

def simulate(
    samples,
    genomes_path,
    reads_path, 
    reference_path=DEFAULT_REFERENCE_PATH,
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

        sample.num_read_pairs = 100
        sample.simulate_reads(simulated.genome_path, reads_path)

        simulated_samples.append(simulated)

    # TODO: pass reads path into function rather than generating inside
    return simulated_samples


def sequence(btb_seq_path, reads_path, results_path):
    # Validate
    if not os.path.exists(btb_seq_path):
        raise Exception("Could not find btb-seq repo at: ", btb_seq_path)

    if not os.path.exists(reads_path):
        raise Exception("Could not find reads path: ", reads_path)

    os.makedirs(results_path, exist_ok=True)

    # Sequence
    btb_seq(btb_seq_path, reads_path, results_path)

    # Result directory
    # TODO: handle when glob does not return a unique path
    return glob.glob(results_path + '/Results_*')[0] + '/'

def benchmark(processed_samples):
    pass

def performance_test(
    output_path,
    btb_seq_path, 
    samples,
    results_path=None,
    simulated_reads_path=None,  
    exist_ok=True, 
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
    
    reads_path = "BLAH"

    simulated_samples = simulate(samples, reads_path)
    processed_samples = sequence(simulated_samples, btb_seq_path, reads_path)
    benchmark(processed_samples)

    return None

    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    output_path = os.path.join(output_path, '')
    if not results_path:
        results_path = output_path
    else:
        results_path = os.path.join(results_path, '')

    # Run simulations if path to simulated reads not provided 
    if not simulated_reads_path:
        simulated_reads_path = simulate(output_path, samples)     
    elif not os.path.isdir(output_path):
        raise Exception("Output path not found")

    # Validate btb_seq_path
    if not os.path.isdir(btb_seq_path):
        raise Exception("Pipeline code repository not found")

    os.makedirs(results_path, exist_ok=exist_ok)
    btb_seq_backup_path = results_path + '/btb-seq/'
    btb_seq_results_path = results_path + '/btb-seq-results/'
    site_stats_path = results_path + '/site-stats/'

    # Paths to mask file and simulated reads directory 
    mask_filepath = btb_seq_backup_path + "references/Mycbovis-2122-97_LT708304.fas.rpt.regions"

    # TODO: handle dwgsim vcf files. Make sure we are taking into account variants it might generate

    # Create Output Directories
    os.makedirs(btb_seq_results_path, exist_ok=exist_ok)
    os.makedirs(site_stats_path, exist_ok=exist_ok)

    # Backup btb-seq code
    # TODO: exclude the work/ subdirectory from this operation.
    #   This could potentially copy large amounts of data
    #   from the work/ directory nextflow generates
    shutil.copytree(btb_seq_path, btb_seq_backup_path, symlinks=True)

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
        samples = quick_samples()
    else:
        samples = standard_samples()

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

class RandomSample2(Sample):
    def __init__(self, 
        reference_path,
            num_snps=16000, 
            num_indels=3898, 
            seed=1, 
            per_base_error_rate="0",
            num_read_pairs = 289994
    ):
        # TODO: Validate
        self.reference_path = reference_path

        self.num_snps = num_snps
        self.num_indels = num_indels
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)
        
    @property
    def name(self):
        return f"{type(self).__name__}-snps{self.num_snps}-indels{self.num_indels}-seed{self.seed}"

    def simulate_genome(self, simulated_genome_path):
        """ Simulated a genome with random SNPs

            TODO: rename simulated_genome_path to simulated_genome_prefix
            Parameters:
                reference_path (str): Path to reference genome
                simulated_genome_path (str): Path to simlated genome
                num_snps (int): Number of random SNPs
                seed (int): Seed value for simulation

            Returns:
                None
        """
        params = [
            "-snp_count", str(self.num_snps),
            "-indel_count", str(self.num_indels),
            "-seed", str(self.seed)
        ]
        self._simulate_genome_base(self.reference_path, simulated_genome_path, params)

        snp_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.SNP.vcf'
        indel_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.INDEL.vcf'
        snp_table_path = simulated_genome_path + '.simulated.refseq2simseq.map.txt'
        genome_path = simulated_genome_path + '.simulated.simseq.genome.fa'

        return simulated_genome.SimulatedGenome(self.name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path)

def benchmark(processed_samples, mask_filepath=DEFAULT_MASK_PATH):
    stats = []
    site_stats = {}
    for sample in processed_samples:
        simulated_snp_path = sample.genome.snp_table_path
        pipeline_snp_path = sample.sequenced.snp_table_path
        pipeline_genome_path = sample.genome.genome_path
        vcf_path = sample.sequenced.vcf_path

        # Performance Stats
        stat = compare_snps.analyse(simulated_snp_path, pipeline_snp_path, pipeline_genome_path, mask_filepath)
        stat["name"] = sample.name
        
        stats.append(stat)

        # Site Statistics at fp/fn/tp positions
        site_stat = compare_snps.site_stats(simulated_snp_path, pipeline_snp_path, vcf_path)

        site_stats[sample.name] = site_stat

    stats_table = pd.DataFrame(stats)
    return stats_table, site_stats

    path = results_path + "stats.csv"
    stats_table.to_csv(path)

# class ProcessedSample:
#     def __init__(self, simulated_sample, sequenced_sample):
#         pass

if __name__ == '__main__':
    # #########
    # results_path = '/home/aaronfishman/temp/vally-1/btb-seq-results/Results_simulated-reads_15Dec21'
    
    

    # print(samples)
    # quit()


    ###### main()
    samples = [RandomSample2('/home/aaronfishman/tinygenome.fas', num_snps=0, num_indels=0)]
        
    genomes_path = '/home/aaronfishman/temp/genomes/'
    reads_path = '/home/aaronfishman/temp/reads/'
    results_path = '/home/aaronfishman/temp/results/'
    btb_seq_path = '/home/aaronfishman/repos/btb-seq/'

    simulated_samples = simulate(samples, genomes_path, reads_path)
    results_path = sequence(btb_seq_path, reads_path, results_path)
    sequenced_samples = sequenced_sample.from_results_dir(results_path)
    processed_samples = processed_sample.from_list(simulated_samples, sequenced_samples)

    stats, site_stats = benchmark(processed_samples)

    print("processed samples", processed_samples)

    print("stats", stats)
    print("site_stats", site_stats)

    a = 1
