from utils import run

class Sample:
    def simulate(self):
        raise NotImplementedError()

class VcfSample(Sample):
    def __init__(self):
        pass

class RandomSample(Sample):
    def __init__(self, num_snps, seed=1):
        self.num_snps = num_snps
        self.seed = seed

    def simulate(self, reference_path, simulated_genome_path):
        simulate_genome_random_snps(reference_path, simulated_genome_path, num_snps=self.num_snps, seed=self.seed)

# reference_path, simulated_genome_path, num_snps=16000, seed=1

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
           "-prefix", simulated_genome_path + "simulated"]
    cmd.extend(params)
    run(cmd)
