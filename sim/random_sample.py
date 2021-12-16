import math

from sample import Sample
from genome import SimulatedGenome
import config

class RandomSample(Sample):
    def __init__(self, 
        reference_path=config.DEFAULT_REFERENCE_PATH,
        num_snps=0,#16000, 
        num_indels=3898, 
        seed=1, 
        per_base_error_rate="0",
        num_read_pairs=config.DEFAULT_NUM_READ_PAIRS
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

        print("ref", self.reference_path)
        print("simulated_genome_path", simulated_genome_path)
        print("params", params)

        self._simulate_genome_base(self.reference_path, simulated_genome_path, params)

        # TODO: DRY with VcfSample
        snp_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.SNP.vcf'
        indel_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.INDEL.vcf'
        snp_table_path = simulated_genome_path + '.simulated.refseq2simseq.map.txt'
        genome_path = simulated_genome_path + '.simulated.simseq.genome.fa'

        return SimulatedGenome(self.name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path)
