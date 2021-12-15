import os

class SimulatedGenome:
    def __init__(self, name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path):
        # Validate
        self.assert_path_exists(genome_path)
        self.assert_path_exists(snp_table_path)
        self.assert_path_exists(snp_vcf_path)
        self.assert_path_exists(indel_vcf_path)

        # Assign
        self.genome_path = genome_path
        self.snp_table_path = snp_table_path
        self.snp_vcf_path = snp_vcf_path
        self.indel_vcf_path = indel_vcf_path
        
        self.name = name

    def assert_path_exists(self, path):
        if not os.path.exists(path):
            raise Exception("Could not find path: ", path)