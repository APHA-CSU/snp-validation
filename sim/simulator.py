from utils import run

class Simulator:
    def __init__(self):
        pass

    def simulate(self, genome_path, read_1_path, read_2_path):
        self.simulate_genome(genome_path)
        self.simulate_reads(genome_path, read_1_path, read_2_path)

    def simulate_reads(self, genome_path, read_1_path, read_2_path):
        raise NotImplementedError("simulate not implemented")

    def simulate_genome(self, genome_path):
        run([
            "simuG.pl",
            "-refseq", genome_path,
            "-snp_count", str(num_snps),
            "-prefix", output_path + "simulated"
        ])

class SuperMutantSimulator(Simulator):
    pass