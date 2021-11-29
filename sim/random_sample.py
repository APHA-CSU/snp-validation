import sample 

class RandomSample(sample.Sample):
    def __init__(self, num_snps=16000, seed=1, per_base_error_rate=sample.PER_BASE_ERROR_RATE):
        self.num_snps = num_snps
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
    
    @property
    def name(self):
        return f"{type(self).__name__}-snps{self.num_snps}-seed{self.seed}"

    def simulate_genome(
        self,
        reference_path, 
        simulated_genome_path,
    ):
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
            "-seed", str(self.seed)
        ]

        self.simug(reference_path, simulated_genome_path, params)
