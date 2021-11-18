from utils import run
import math 

import os

class Sample:
    per_base_error_rate = "0"
    num_read_pairs = 150000

    def simulate_genome(self):
        raise NotImplementedError()

    @property
    def name(self):
        # TODO: this is ugly, make this human-readable
        return f"unnamed-{str(id(self))}"

    def simulate_reads(self, simulated_genome_path, simulated_reads_path):
        # # TODO: explicitly path fasta path to simulate
        fasta_path = simulated_genome_path + self.name + '.simulated.simseq.genome.fa'
        simulate_reads(
            fasta_path, 
            simulated_reads_path, 
            sample_name=self.name, 
            per_base_error_rate=self.per_base_error_rate,
            num_read_pairs=self.num_read_pairs
        )

class NamedSample(Sample):
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name


class VcfSample(Sample):
    def __init__(self, 
        predef_snp_path, 
        seed=1, 
        per_base_error_rate="0",
        num_read_pairs = 10,
    ):
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)

        # TODO: check predef_snp_path exists
        self.predef_snp_path = predef_snp_path

    @property
    def name(self):
        # extract filename
        #TODO: this will fail if the extension is not 3 charachters long
        basename = os.path.basename(self.predef_snp_path)[:-4]
        return f"{type(self).__name__}-{basename}-seed{self.seed}-error{self.per_base_error_rate}-readpairs-{self.num_read_pairs}"

    def simulate_genome(self, reference_path, simulated_genome_path):
        simulate_genome_from_vcf(reference_path, simulated_genome_path, self.predef_snp_path, seed=1)

class RandomSample(Sample):
    def __init__(self, num_snps=16000, seed=1, per_base_error_rate="0"):
        self.num_snps = num_snps
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
    
    @property
    def name(self):
        return f"{type(self).__name__}-snps{self.num_snps}-seed{self.seed}"

    def simulate_genome(self, reference_path, simulated_genome_path):
        simulate_genome_random_snps(reference_path, simulated_genome_path, num_snps=self.num_snps, seed=self.seed)

# TODO: move these functions inside of the classes

def simulate_genome_random_snps(reference_path, simulated_genome_path, num_snps=16000, seed=1):
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
    params = ["-snp_count", str(num_snps),
                "-seed", str(seed)]
    simulate_genome(reference_path, simulated_genome_path, params)

def simulate_genome_from_vcf(reference_path, simulated_genome_path, predef_snp_path, seed=1):
    """ Simulated a genome with random SNPs

        TODO: rename simulated_genome_path to simulated_genome_prefix

        Parameters:
            reference_path (str): Path to reference genome
            simulated_genome_path (str): Path to simlated genome
            predef_snps_path (str): Path to snippy generated VCF file.
            seed (int): Seed value for simulation

        Returns:
            None
    """
    params = ["-snp_vcf", predef_snp_path, 
              #"-indel_vcf", predef_snp_path # tells simuG to also simulate predefined indels.
              "-seed", str(seed)]
    simulate_genome(reference_path, simulated_genome_path, params)

def simulate_genome(reference_path, simulated_genome_path, params):
    cmd = ["simuG.pl",
           "-refseq", reference_path,
           "-prefix", simulated_genome_path + ".simulated"]
    cmd.extend(params)
    run(cmd)

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
    per_base_error_rate="0.001-0.01" # TODO: default to Ele's reccomendation? 0.001-0.01
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